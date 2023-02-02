import networkx as nx
import math
import gurobipy as gp
from gurobipy import GRB

import matplotlib.pyplot as plt
from pyproj import Proj

from util import update_attributes, get_k_L_U, get_epsg
from separation import labeling_contiguity_callback

def sort_by_second(val):
    return val[1]

def sq_eucl_dist(x1,y1,x2,y2):
    return (x1-x2)**2 + (y1-y2)**2

def find_ordering(DG, source, level):
    if level == 'tract':
        n = DG.number_of_nodes()
        for (i,j) in DG.edges:
            gi = DG.nodes[i]['GEOID20'][0:5]
            gj = DG.nodes[j]['GEOID20'][0:5]
            if gi == gj:  # endpoints belong to same county
                DG[i][j]['weight'] = 1
            else:         # endpoints belong to different counties
                DG[i][j]['weight'] = n
    elif level == 'block':
        n = DG.number_of_nodes()
        n_squared = n * n
        for (i,j) in DG.edges:
            gi = DG.nodes[i]['GEOID20']
            gj = DG.nodes[j]['GEOID20']
            if gi[0:11] == gj[0:11]: # endpoints belong to same tract
                DG[i][j]['weight'] = 1
            elif gi[0:5] == gj[0:5]:  # endpoints belong to different tracts, but same county
                DG[i][j]['weight'] = n 
            else: # endpoints belong to different counties
                DG[i][j]['weight'] = n_squared
    else:
        for (i,j) in DG.edges:
            DG[i][j]['weight'] = 1
            
    dist = nx.shortest_path_length(DG,source=source,weight='weight')
    vertex_dist_pairs = [ (i,dist[i]) for i in dist.keys() ]
    vertex_dist_pairs.sort(key=sort_by_second)
    return [ v for (v,dist) in vertex_dist_pairs ]


def find_positions(G, ordering):
    pos = { v : None for v in G.nodes }
    for p in range(len(ordering)):
        v = ordering[p]
        pos[v] = p
    return pos


def detail(G, L, U, k, county_geoid_support=None, time_limit=600):
    
    assert k >= 1
    if k == 1:
        assert nx.is_connected( G )
        total_population = sum( G.nodes[i]['TOTPOP'] for i in G.nodes )
        assert total_population >= L
        assert total_population <= U
        if county_geoid_support is not None:
            for i in G.nodes:
                g = G.nodes[i]['GEOID20'][0:5]
                assert 0 in county_geoid_support[g] 
        return [ list(G.nodes) ]
    
    # tract- or block-level?
    v = nx.utils.arbitrary_element(G.nodes)
    if len( G.nodes[v]['GEOID20'] ) == 15:
        print("Trying block-level instance.")
        level = 'block'
    else:
        print("Trying tract-level instance.")
        level = 'tract'
    
    # build model
    m = gp.Model()
    m.Params.OutputFlag = 0
    x = m.addVars( G.nodes, k, vtype=GRB.CONTINUOUS )
    
    if county_geoid_support is not None:
        
        # initially disallow all assignments
        for i in G.nodes:
            for j in range(k):
                x[i,j].UB = 0

        # permit assignments in the support
        for i in G.nodes:
            g = G.nodes[i]['GEOID20'][0:5]
            for j in county_geoid_support[g]:
                x[i,j].UB = 1
                
    # assignment constraints
    m.addConstrs( gp.quicksum( x[i,j] for j in range(k) ) == 1 for i in G.nodes )
    
    # population balance
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in G.nodes ) >= L for j in range(k) )
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in G.nodes ) <= U for j in range(k) )
    
    max_iterations = 100
    print("Iteration \t Objective \t iter_time")
    
    for iteration in range(max_iterations):
    
        if iteration == 0:
            old_total_cost = -1               # dummy value
            m.setObjective( 0, GRB.MINIMIZE ) # dummy objective
        else:
            old_total_cost = total_cost
            m.setObjective( gp.quicksum( G.nodes[i]['TOTPOP'] * sq_eucl_dist(G.nodes[i]['X'],G.nodes[i]['Y'],means_x[j],means_y[j]) * x[i,j] for i in G.nodes for j in range(k) ), GRB.MINIMIZE )

        m.optimize()
        print(iteration,"\t",'{0:.2f}'.format(m.objVal),"\t",'{0:.2f}'.format(m.runtime))
    
        total_cost = m.objVal
        clusters = [ [ i for i in G.nodes if x[i,j].x > 0.5 ] for j in range(k) ]

        if total_cost == old_total_cost:
            break
        else:
            # get next means
            population = [ sum( G.nodes[i]['TOTPOP'] * x[i,j].x for i in G.nodes ) for j in range(k) ]
            means_x = [ sum( G.nodes[i]['TOTPOP'] * G.nodes[i]['X'] * x[i,j].x for i in G.nodes ) / population[j] for j in range(k) ]
            means_y = [ sum( G.nodes[i]['TOTPOP'] * G.nodes[i]['Y'] * x[i,j].x for i in G.nodes ) / population[j] for j in range(k) ]
    
    # convert to binary vars
    for i in G.nodes:
        for j in range(k):
            if x[i,j].UB > 0.5:
                x[i,j].vtype = GRB.BINARY
    
    # parameters
    m.Params.TimeLimit = time_limit
    m.Params.MIPGap = 0.10 
    m.Params.FeasibilityTol = 1e-7
    m.Params.IntFeasTol = 1e-7
    m.Params.OutputFlag = 1
    
    ########################
    # Add DAG constraints
    ########################
    
    # first, identify roots
    roots = list()
    for j in range(k):
        min_dist = None
        min_dist_node = None
        for i in G.nodes:
            g = G.nodes[i]['GEOID20'][0:5]
            if j not in county_geoid_support[g]:
                continue
            dist = sq_eucl_dist( G.nodes[i]['X'], G.nodes[i]['Y'], means_x[j], means_y[j] )
            if min_dist is None or dist < min_dist:
                min_dist = dist
                min_dist_node = i
        roots.append(min_dist_node)
    
    # now add constraints
    DG = nx.DiGraph(G)
    c = [ None for j in range(k) ]
    for j in range(k):
        root = roots[j]
        x[root,j].LB = 1
        order = find_ordering(DG, root, level)
        pos = find_positions(G, order)
        c[j] = m.addConstrs( x[i,j] <= gp.quicksum( x[p,j] for p in G.neighbors(i) if pos[p] < pos[i] ) for i in G.nodes if i != root )
    
    print("Trying DAG model, with roots =",roots)
    m.optimize()
    
    if m.solCount > 0:
        districts = [ [ i for i in G.nodes if x[i,j].x > 0.5 ] for j in range(k) ]
        return districts

    m.remove(c)
    for j in range(k):
        root = roots[j]
        x[root,j].LB = 0
        
    ##################################################################
    # if working with a block-level graph, 
    # then permit at most k-1 tracts to be split.
    # moreover, let its blocks only be assigned to two districts
    ##################################################################
    
    if level == 'block':
        tract_geoids = { G.nodes[i]['GEOID20'][0:11] for i in G.nodes }
        
        # ttd[g,j]=1 if some of tract with geoid g is assigned to district j
        ttd = m.addVars( tract_geoids, k, vtype=GRB.BINARY )
        
        # only permit at most # tracts + (k-1) assignments
        m.addConstr( gp.quicksum(ttd) <= len(tract_geoids) + k - 1 )
        
        # split each tract at most once
        m.addConstrs( gp.quicksum( ttd[g,j] for j in range(k) ) <= 2 for g in tract_geoids )
        
        # block i can only be assigned to district j if its tract is 
        m.addConstrs( x[i,j] <= ttd[G.nodes[i]['GEOID20'][0:11],j] for i in G.nodes for j in range(k) )
    
    #####################################
    # Add a,b-separator constraints
    #####################################
    
    print("Trying a,b-separator model.")
    m.Params.LazyConstraints = 1
    m._callback = labeling_contiguity_callback
    m._DG = DG
    m._x = x
    m._parts = k
    m.optimize(m._callback)

    if m.solCount > 0:
        districts = [ [ i for i in G.nodes if m._x[i,j].x > 0.5 ] for j in range(k) ]
    else:
        districts = None
        
    return districts