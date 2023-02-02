import networkx as nx
import math
import gurobipy as gp
from gurobipy import GRB  

from separation import labeling_contiguity_callback

def sketch(G, L, U, k):
    
    m = gp.Model()
    m.Params.OutputFlag = 0
    
    # x[i,j]=1 if county i is assigned to district j (in whole or in part)
    x = m.addVars( G.nodes, k, vtype=GRB.BINARY )
    
    # y[e,j]=1 if edge e is assigned to district j (in whole or in part)
    y = m.addVars( G.edges, k, vtype=GRB.BINARY )
    
    # ysum[e]=1 if edge e is used at all
    ysum = m.addVars( G.edges, vtype=GRB.BINARY )
    m.addConstrs( ysum[u,v] == gp.quicksum( y[u,v,j] for j in range(k) ) for u,v in G.edges )
    
    # z[i,j] = the proportion of county i that is assigned to district j
    z = m.addVars( G.nodes, k )
    
    # s[i] = the number of splits permitted in county i
    s = m.addVars( G.nodes, vtype=GRB.INTEGER )
    for i in G.nodes:
        s[i].LB = math.ceil( G.nodes[i]['TOTPOP'] / U ) - 1  
           
    # indegree extended formulation
    m.addConstrs( gp.quicksum( x[i,j] for i in G.nodes ) - gp.quicksum( y[u,v,j] for u,v in G.edges ) <= 1 for j in range(k) )
    
    # if pick edge, then must pick endpoints
    for u,v in G.edges:
        m.addConstrs( y[u,v,j] <= x[u,j] for j in range(k) )
        m.addConstrs( y[u,v,j] <= x[v,j] for j in range(k) )
        
    # pick edge at most once
    m.addConstrs( gp.quicksum( y[u,v,j] for j in range(k) ) <= 1 for u,v in G.edges )
    
    # if pick endpoints, then must pick edge
    for u,v in G.edges:
        m.addConstrs( y[u,v,j] >= x[u,j] + x[v,j] - 1 for j in range(k) )

    # assignment constraints
    m.addConstrs( gp.quicksum( z[i,j] for j in range(k) ) == 1 for i in G.nodes )
        
    # population balance
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * z[i,j] for i in G.nodes ) >= L for j in range(k) )
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * z[i,j] for i in G.nodes ) <= U for j in range(k) )

    # coupling constraints
    m.addConstrs( z[i,j] <= x[i,j] for i in G.nodes for j in range(k) )

    # splitting constraints
    m.addConstrs( gp.quicksum( x[i,j] for j in range(k) ) == s[i] + 1 for i in G.nodes )
    m.addConstr( gp.quicksum(s) == k-1 )
        
    # eliminate some symmetry
    max_pop = max( G.nodes[i]['TOTPOP'] for i in G.nodes )
    max_pop_node = [ i for i in G.nodes if G.nodes[i]['TOTPOP'] == max_pop ][0]
    m.addConstrs( z[max_pop_node,j] >= z[max_pop_node,j+1] for j in range(k-1) )
    lb = math.ceil( max_pop / U )
    for j in range(lb):
        x[max_pop_node,j].LB = 1
        
    # contiguity constraints
    m.Params.LazyConstraints = 1
    m._DG = nx.DiGraph(G)
    m._parts = k
    m._x = x
    
    # first, maximize preserved edges
    m.setObjective( gp.quicksum(y), GRB.MAXIMIZE )
    #m.setObjective( gp.quicksum( G.edges[u,v]['shared_perim'] * y[u,v,j] for u,v in G.edges for j in range(k) ), GRB.MAXIMIZE )
    
    m._callback = labeling_contiguity_callback
    m.Params.MIPGap = 0.10
    m.optimize(m._callback)
    
    if m.solCount > 0:
    
        # fix the number of preserved edges
        obj = round( m.objVal )
        m.addConstr( gp.quicksum(y) >= obj )
        
        #obj = m.objVal - 0.0001
        #m.addConstr( gp.quicksum( G.edges[u,v]['shared_perim'] * y[u,v,j] for u,v in G.edges for j in range(k) ) >= obj )

        # c[j,t] = 1 if district j touches t counties
        c = m.addVars( k, range(1,G.number_of_nodes()+1), vtype=GRB.BINARY )
        m.addConstrs( gp.quicksum( c[j,t] for t in range(1,G.number_of_nodes()+1) ) == 1 for j in range(k) )
        m.addConstrs( gp.quicksum( t * c[j,t] for t in range(1,G.number_of_nodes()+1) ) == gp.quicksum( x[i,j] for i in G.nodes ) for j in range(k) )

        # minimize sum of squares of counties touched 
        m.setObjective( gp.quicksum( t * t * c[j,t] for j in range(k) for t in range(1,G.number_of_nodes()+1) ), GRB.MINIMIZE )
        m._callback = labeling_contiguity_callback
        m.optimize(m._callback)

        xsoln = { i : { j : x[i,j].x for j in range(k) if x[i,j].x > 0 } for i in G.nodes }
        ysoln = { (u,v) : { j : y[u,v,j].x for j in range(k) if y[u,v,j].x > 0 } for (u,v) in G.edges }
        zsoln = { i : { j : z[i,j].x for j in range(k) if z[i,j].x > 0 } for i in G.nodes }
    else:
        xsoln = None
        ysoln = None
        zsoln = None
        
    return (xsoln, ysoln, zsoln)
