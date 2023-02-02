import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from itertools import combinations
import math
import random
import time

from separation import labeling_contiguity_callback

def sort_by_second(val):
    return val[1]

def find_ordering(G):
    V_with_population = [ (i, G.nodes[i]['TOTPOP']) for i in G.nodes ]
    V_with_population.sort( key=sort_by_second, reverse=True )
    return [ v for (v,p) in V_with_population ]

def construct_position(ordering):
    position = { i : -1 for i in ordering }
    for p in range(len(ordering)):
        v = ordering[p]
        position[v] = p
    return position

# key:    ( tuple(sorted(G.nodes)), k, num_clusters )
# value:  (clusters, sizes, cut_edge_LB, cut_edge_UB)
min_cut_cache = dict()
        
def min_cut_county_clustering_via_labeling(G, L, U, k, num_clusters, time_limit=3600, carve=False, initial_clusters=None, initial_sizes=None, randomized_carving=False, cache=True, verbose=False):
    
    # check if we've already solved this instance before.
    # if so, return our cached solution!
    key = ( tuple(sorted(G.nodes)), k, num_clusters )
    if cache and key in min_cut_cache:
        #print("Pulling min_cut solution from cache.")
        return min_cut_cache[key]
    
    # when carving, we only seek 2 clusters!
    assert num_clusters == 2 or carve == False
    
    DG = nx.DiGraph(G) # bidirected version of G
    m = gp.Model()
    if not verbose:
        m.Params.OutputFlag = 0
    m.Params.TimeLimit = time_limit

    # VARIABLES
    # x[i,j]=1 if vertex i is assigned to cluster j
    x = m.addVars(DG.nodes, num_clusters, vtype=GRB.BINARY) 

    # y[j] = # of districts in cluster j
    y = m.addVars(num_clusters, vtype=GRB.INTEGER, lb=1)
    
    # z[u,v,j] = 1 if directed edge (u,v) is cut because u->j but v!->j
    z = m.addVars(DG.edges, num_clusters, vtype=GRB.BINARY)
    
    # is_cut[u,v] = 1 if undirected edge {u,v} is cut
    is_cut = m.addVars(G.edges, vtype=GRB.BINARY)
    
    # BASE CONSTRAINTS
    # Each county i assigned to one cluster
    m.addConstrs( gp.quicksum( x[i,j] for j in range(num_clusters) ) == 1 for i in DG.nodes )

    # Cluster "sizes" should sum to k
    m.addConstr( gp.quicksum( y ) == k )

    # Population balance: population of cluster j should be an integer multiple of [L,U] 
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in DG.nodes ) <= U * y[j] for j in range(num_clusters) )
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in DG.nodes ) >= L * y[j] for j in range(num_clusters) )
    
    # CUT EDGE CONSTRAINTS
    m.addConstrs( x[u,j] - x[v,j] <= z[u,v,j] for u,v in DG.edges for j in range(num_clusters) )
    m.addConstrs( is_cut[u,v] == gp.quicksum( z[u,v,j] for j in range(num_clusters) ) for u,v in G.edges )
    m.addConstrs( is_cut[u,v] == gp.quicksum( z[v,u,j] for j in range(num_clusters) ) for u,v in G.edges )
                
    # SYMMETRY BREAKING, via diagonal fixing
    # Find vertex ordering
    ordering = find_ordering(G)

    for p in range(G.number_of_nodes()):
        i = ordering[p]
        for j in range(p+1,num_clusters):
            x[i,j].UB = 0
            
    # MIP START
    if initial_clusters is not None:

        assert len(initial_clusters) == num_clusters
        assert len(initial_sizes) == num_clusters
        
        root_found = [ False for j in range(num_clusters) ]
        labeling = { i : j for j in range(num_clusters) for i in initial_clusters[j] } 
        p = 0
        for i in ordering:
            j = labeling[i]
            if not root_found[j]:
                root_found[j] = True
                y[p].start = initial_sizes[j]
                for v in initial_clusters[j]:
                    x[v,p].start = 1
                p += 1

    # CONTIGUITY CONSTRAINTS
    m.Params.LazyConstraints = 1
    m._DG = DG
    m._parts = num_clusters
    m._x = x
    m._callback = labeling_contiguity_callback
    
    if carve:
        # *temporarily* use large MIP gap, just to determine feasibility
        m.Params.MIPGap = 1.00
    
    # MINIMIZE CUT EDGES
    m.setObjective( gp.quicksum(is_cut), GRB.MINIMIZE )
    m.Params.IntFeasTol = 1.e-7
    m.Params.FeasibilityTol = 1.e-7
    m.optimize( m._callback )
    
    # NO SOLUTION FOUND? EXIT.
    if m.solCount == 0:
        clusters = [ list(G.nodes) ]
        sizes = [ k ]
        cut_edge_LB = 0
        cut_edge_UB = G.number_of_edges() + 1
        
        # add solution to our cache and return
        min_cut_cache[key] = (clusters, sizes, cut_edge_LB, cut_edge_UB)
        return min_cut_cache[key]
    
    if carve:
        m.Params.MIPGap = 0.00
        # give the complement solution as a warm start
        for j in range(2):
            y[j].start = k - round( y[j].x )
            for i in G.nodes:
                m._x[i,j].start = 1 - round( m._x[i,j].x )
        
        # undo diagonal fixing
        for p in range(G.number_of_nodes()):
            i = ordering[p]
            for j in range(p+1,num_clusters):
                x[i,j].UB = 1
        
        # try each possible size for first cluster, and return
        # the smallest feasible one, with fewest cut edges
        if randomized_carving:
            for i,j in G.edges:
                is_cut[i,j].obj = random.random() # random number between 0 and 1
        else:
            m.setObjective( gp.quicksum(is_cut), GRB.MINIMIZE )
        for size in range(1,k-1):
            y[0].LB = size
            y[0].UB = size
            m._callback = labeling_contiguity_callback
            m.optimize( m._callback )
            if m.solCount > 0:
                break
        
    # Retrieve solution
    clusters = [ [ i for i in G.nodes if m._x[i,j].x > 0.5 ] for j in range(num_clusters) ]
    sizes = [ round( y[j].x ) for j in range(num_clusters) ]
    cut_edge_LB = round( m.objBound )
    cut_edge_UB = round( m.objVal )

    # add solution to our cache and return
    if cache:
        min_cut_cache[key] = (clusters, sizes, cut_edge_LB, cut_edge_UB)
    return (clusters, sizes, cut_edge_LB, cut_edge_UB)

# Is there an integer size for which 
#   L*size <= population <= U*size ?
#
def is_within_window(population, L, U):
    min_size = math.ceil( population / U )
    max_size = math.floor( population / L )
    return min_size <= max_size

# Is it possible to build a county cluster 
#  that contains all nodes from 'required_nodes'?
#
def is_cluster_possible(G, required_nodes, L, U, k, verify=True):
    (cluster, size) = find_cluster(G, required_nodes, L, U, k, verify) 
    return cluster is not None
    
def find_cluster(G, required_nodes, L, U, k, verify=True):
    
    # To make contiguity easy, this function assumes that
    #   required_nodes is a connected dominating set.
    if verify:
        assert nx.is_connected( G.subgraph( required_nodes ) )
        assert nx.is_dominating_set( G, required_nodes )
    
    m = gp.Model()
    m.Params.OutputFlag = 0
    x = m.addVars( G.nodes, vtype=GRB.BINARY )
    for j in required_nodes:
        x[j].LB = 1
    
    z = m.addVar( vtype=GRB.INTEGER )
    z.LB = 1
    z.UB = k
    
    m.addConstr( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i] for i in G.nodes ) >= L * z )
    m.addConstr( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i] for i in G.nodes ) <= U * z )
    
    m.optimize()
    if m.solCount > 0:
        cluster = [ i for i in G.nodes if x[i].x > 0.5 ]
        size = round(z.x)
        return ( cluster, size ) 
    else:
        return (None, None)


# key:    ( tuple(sorted(G.nodes)), k)
# value:  (clusters, sizes, cluster_UB) 
max_cluster_cache = dict()

def max_county_clustering_via_hess(G, L, U, k, time_limit=3600, initial_clusters=None, initial_sizes=None, valid_inequalities=True, cache=True, verbose=False): 
    
    # use cluster_UB to exit early?
    cluster_UB = k - sum( math.floor( G.nodes[i]['TOTPOP'] / (U+1) ) for i in G.nodes )
    if initial_sizes is not None and len(initial_sizes) == cluster_UB:
        return (initial_clusters, initial_sizes, cluster_UB)
    
    # use cache to exit early?
    key = ( tuple(sorted(G.nodes)), k)
    if cache and key in max_cluster_cache:
        #print("Pulling max_cluster solution from cache.")
        return max_cluster_cache[key]
    
    # setup    
    DG = nx.DiGraph(G) # bidirected version of G
    m = gp.Model()
    if not verbose:
        m.Params.OutputFlag = 0
    m.Params.TimeLimit = time_limit

    # x[i,j]=1 if vertex i is assigned to (cluster centered at) vertex j
    x = m.addVars(DG.nodes, DG.nodes, vtype=GRB.BINARY) 

    # y[j] = # of districts in cluster centered at vertex j
    y = m.addVars(DG.nodes, vtype=GRB.INTEGER)
    
    # BASE CONSTRAINTS
    # Each county i assigned to one cluster
    m.addConstrs( gp.quicksum( x[i,j] for j in DG.nodes ) == 1 for i in DG.nodes )

    # Cluster "sizes" should sum to k
    m.addConstr( gp.quicksum( y ) == k )

    # Population balance: population of cluster j should be an integer multiple of [L,U] 
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in DG.nodes ) <= U * y[j] for j in DG.nodes )
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in DG.nodes ) >= L * y[j] for j in DG.nodes )

    # Enforce relationships between x[i,j] and x[j,j]
    m.addConstrs( x[i,j] <= x[j,j] for i in DG.nodes for j in DG.nodes if i != j )
    
    # CONTIGUITY CONSTRAINTS
    # f[j,u,v] tells how much flow (from source j) is sent across arc (u,v)
    f = m.addVars(DG.nodes, DG.edges)
    M = DG.number_of_nodes() - 1
    m.addConstrs( gp.quicksum( f[j,u,j] for u in DG.neighbors(j) ) == 0 for j in DG.nodes )
    m.addConstrs( gp.quicksum( f[j,u,i] - f[j,i,u] for u in DG.neighbors(i) ) == x[i,j] for i in DG.nodes for j in DG.nodes if i != j )
    m.addConstrs( gp.quicksum( f[j,u,i] for u in DG.neighbors(i) ) <= M * x[i,j] for i in DG.nodes for j in DG.nodes if i != j )
    
    # Find vertex ordering
    ordering = find_ordering(G)
    position = construct_position( ordering )

    # Diagonal fixing
    for i in G.nodes:
        for j in G.nodes:
            if position[i] < position[j]:
                x[i,j].UB = 0
    m.update()
    
    # Add Chvatal-Gomory-like valid inequalities (t=1)
    if valid_inequalities:
        for j in G.nodes:
            if y[j].UB > 0.5:
                m.addConstr( gp.quicksum( math.floor( G.nodes[i]['TOTPOP'] / (U+1) ) * x[i,j] for i in G.nodes ) <= y[j] - x[j,j] )
    
    # Add cluster-separator inequalities (special cases, anyway)
    # for each vertex b, try a=b and S=N(b)
    if valid_inequalities:
        for b in G.nodes:
            if y[b].UB > 0.5 and not is_within_window(G.nodes[b]['TOTPOP'], L, U):
                m.addConstr( x[b,b] <= sum( x[i,b] for i in G.neighbors(b) if position[i] > position[b] ) )
    
    # for each edge {a,b}, try S=N({a,b})
    if valid_inequalities:
        for a,b in DG.edges:
            if position[a] < position[b] or y[b].UB < 0.5: # a cannot be assigned to b anyway..
                continue
            if not is_within_window(G.nodes[a]['TOTPOP']+G.nodes[b]['TOTPOP'], L, U):
                S = nx.node_boundary(G, [a,b])
                S = [ i for i in S if position[i] > position[b] ]
                m.addConstr( x[a,b] <= sum( x[i,b] for i in S ) )
    
    # for each vertex b, try a=b and all S=N(component) for all components "around" b. 
    if valid_inequalities:
        for b in G.nodes:
            if y[b].UB < 0.5:
                continue
            for size in range(2,1+G.degree(b)):
                for comb in combinations(G.neighbors(b),size):
                    component = list(comb) + [b]
                    component = [ i for i in component if position[i] >= position[b] ]
                    if not is_cluster_possible( G.subgraph(component), [b], L, U, k, verify=False ):
                        S = nx.node_boundary(G, component)
                        S = [ i for i in S if position[i] > position[b] ]
                        m.addConstr( x[b,b] <= sum( x[i,b] for i in S ) )
                    
    # Provide MIP start
    if initial_clusters is not None:
        assert len(initial_clusters) == len(initial_sizes)
        for p in range(len(initial_clusters)):
            cluster = initial_clusters[p]
            min_pos = min( position[i] for i in cluster )
            j = [ i for i in cluster if position[i] == min_pos ][0]
            y[j].start = initial_sizes[p]
            for i in cluster:
                x[i,j].start = 1 
                
    # Parameter settings
    m.Params.IntFeasTol = 1.e-7
    m.Params.FeasibilityTol = 1.e-7
    for j in DG.nodes:
        x[j,j].BranchPriority = 1
        y[j].BranchPriority = 1
    
    # Maximize # of clusters
    m.setObjective( gp.quicksum( x[j,j] for j in G.nodes ), GRB.MAXIMIZE )
    m.optimize()
    cluster_UB = round( m.objBound )
    
    # Retrieve solution
    if m.solCount > 0:
        centers = [ j for j in DG.nodes if x[j,j].x > 0.5 ]
        clusters = [ [ i for i in G.nodes if x[i,j].x > 0.5 ] for j in centers ]
        sizes = [ round( y[j].x ) for j in centers ]
        # double check the populations and sizes
        assert sum( sizes ) == k
        for p in range(len(sizes)):
            pop = sum( G.nodes[i]['TOTPOP'] for i in clusters[p] )
            assert pop >= L * sizes[p] and pop <= U * sizes[p]
    else:
        clusters = [ list(G.nodes) ]
        sizes = [k]
        
    # add solution to our cache and return
    if cache:
        max_cluster_cache[key] = (clusters, sizes, cluster_UB)  
    return (clusters, sizes, cluster_UB)  


def recom_t_opt(G, clusters, sizes, L, U, t, try_max_cluster=True, try_min_cut=True):
    
    assert try_max_cluster or try_min_cut # otherwise, why are you here???
    
    num_clusters = len( clusters )
    for comb in combinations( range(num_clusters), t ):
        
        improved_clusters = False
        improved_cut_edges = False
        
        # if the cluster union is disconnected, then skip it.
        cluster_union = [ i for p in comb for i in clusters[p] ]
        if not nx.is_connected( G.subgraph(cluster_union) ):
            continue
        
        # incumbents
        temp_clusters = [ clusters[p] for p in comb ]
        temp_sizes = [ sizes[p] for p in comb ]
        cluster_size = sum( temp_sizes )
        num_clusters = t
        
        ###############################
        # increase number of clusters?
        ###############################
        if try_max_cluster:
            
            (hess_clusters, hess_sizes, cluster_UB) = max_county_clustering_via_hess(G.subgraph(cluster_union), L, U, cluster_size, initial_clusters=temp_clusters, initial_sizes=temp_sizes) 

            num_clusters = len(hess_clusters)

            if num_clusters > t:
                cut_cost = number_of_cut_edges(G.subgraph(cluster_union),hess_clusters)
                cut_cost -= number_of_cut_edges(G.subgraph(cluster_union),temp_clusters)
                print("clusters +=", num_clusters-t,"( w/ cut edges +=",cut_cost,")")
                improved_clusters = True
                temp_clusters = hess_clusters
                temp_sizes = hess_sizes

        ################################
        # decrease number of cut edges?
        ################################
        if try_min_cut:
            
            old_cut_edges = number_of_cut_edges( G.subgraph(cluster_union), temp_clusters )

            (temp_clusters, temp_sizes, cut_edge_LB, cut_edge_UB) = min_cut_county_clustering_via_labeling( G.subgraph(cluster_union), L, U, cluster_size, num_clusters, time_limit=60, initial_clusters=temp_clusters, initial_sizes=temp_sizes)

            if cut_edge_UB < old_cut_edges:
                print("cut edges -=",old_cut_edges-cut_edge_UB)
                improved_cut_edges = True
        
        ##############################################
        # found a local search move? if so, update clusters and recurse
        ##############################################
        if improved_clusters or improved_cut_edges:
            
            # write over the existing clusters
            for p in range(t):
                j = comb[p]
                clusters[j] = temp_clusters[p]
                sizes[j] = temp_sizes[p]
            
            # append the 'extra' clusters that this soln found (if any)
            if improved_clusters:
                for p in range(t, num_clusters ):
                    clusters.append( temp_clusters[p] )
                    sizes.append( temp_sizes[p] )
            
            return recom_t_opt(G, clusters, sizes, L, U, t, try_max_cluster, try_min_cut)
    
    # Failed to find an improving move; return the input clusters/sizes.
    return ( clusters, sizes )


def number_of_cut_edges( G, clusters ):
    labeling = { i : j for j in range(len(clusters)) for i in clusters[j] }
    return sum( 1 for i,j in G.edges if labeling[i] != labeling[j] )


def carve_cluster(G, L, U, k, randomized_carving=False):
    (clusters, sizes, cut_edge_LB, cut_edge_UB) = min_cut_county_clustering_via_labeling(G, L, U, k, num_clusters=2, carve=True, randomized_carving=randomized_carving, cache=not randomized_carving)
    return (clusters[0], sizes[0])
    

def carve_heuristic(G, L, U, k, randomized_carving=False):
    
    clusters = list()
    sizes = list()
    R = set(G.nodes)   # remaining vertices
    kr = k              # remaining size
    
    while kr > 0:
        (cluster, size) = carve_cluster( G.subgraph(R), L, U, kr, randomized_carving )
        if len(clusters)==0:
            print("carved cluster sizes = ",end="")
            
        clusters.append(cluster)
        sizes.append(size)
        R -= set( cluster )
        kr -= size
        print(size,end=', ')
        
    print("\ncarved LB =", len(clusters))
    print("carved cut edges =",number_of_cut_edges(G, clusters))
    return (clusters, sizes)

# converts our double objective into a single value
#   1st objective: max # clusters
#   2nd objective: min # cut edges
def objective(G, clusters):
    M = G.number_of_edges() + 1
    return M * len(clusters) - number_of_cut_edges(G, clusters)

def max_cluster_main(G, L, U, k, max_t=4, restarts=3, time_limit=3600):   
    
    results = dict()
    
    # incumbents
    (clusters, sizes) = ( [ list(G.nodes) ], [k] )
    cluster_UB = k - sum( math.floor( G.nodes[i]['TOTPOP'] / (U+1) ) for i in G.nodes )
    results['initial_UB'] = cluster_UB
    print("Initially, cluster_UB =",cluster_UB)
    
    if max_t > cluster_UB - 1:
        print("No need for t-opt local search, with t =",max_t,"; reducing to",cluster_UB - 1)
        max_t = cluster_UB - 1
        
    # apply carving construction heuristic and t-opt recom local search
    start_time = time.time()
    for iteration in range(restarts):
        print("\n****************************")
        print("Heuristic iteration #",iteration)
        print("****************************")
        (iclusters, isizes) = carve_heuristic(G, L, U, k, randomized_carving=True)
        for t in range(2, max_t+1):
            (iclusters, isizes) = recom_t_opt(G, iclusters, isizes, L, U, t)
            print("t =",t,"-> #clusters, #cut edges =",len(iclusters),number_of_cut_edges(G, iclusters) )
        if objective(G, iclusters) > objective(G, clusters):
            (clusters, sizes) = (iclusters, isizes)
            print("new incumbent!")

    results['heuristic_time'] = '{0:.2f}'.format( time.time() - start_time )
    results['heuristic_num_clusters'] = len(clusters)
    results['heuristic_num_cut_edges'] = number_of_cut_edges(G, clusters)
    results['heuristic_iterations'] = restarts
    
    print("\n********************************************************")
    print("After local search, # clusters, #cut edges =",results['heuristic_num_clusters'],results['heuristic_num_cut_edges'] )
    print("********************************************************\n")
    
    # warm start our max_cluster MIP 
    start_time = time.time()
    (hess_clusters, hess_sizes, cluster_UB) = max_county_clustering_via_hess(G, L, U, k, time_limit=time_limit, initial_clusters=clusters, initial_sizes=sizes, valid_inequalities=True, cache=False, verbose=True)
    results['MIP_time'] = '{0:.2f}'.format( time.time() - start_time )
    
    print("********************************************************")
    print("MIP gives #clusters, #cut edges =",len(hess_clusters),number_of_cut_edges(G, hess_clusters) )
    print("********************************************************\n")

    # last, make MIP-generated clusters more compact (only if they are different from before)
    start_time = time.time()
    if len(hess_clusters) > results['heuristic_num_clusters']:
        (clusters, sizes) = (hess_clusters, hess_sizes)
        for t in range(2, max_t+1):
            (clusters, sizes) = recom_t_opt(G, clusters, sizes, L, U, t, try_max_cluster=False)
            print("t =",t,"-> #clusters, #cut edges =",len(clusters),number_of_cut_edges(G, clusters) )
            
    results['cleanup_time'] = '{0:.2f}'.format( time.time() - start_time )

    results['clusters'] = clusters
    results['sizes'] = sizes
    results['num_clusters'] = len(clusters)
    results['num_cut_edges'] = number_of_cut_edges(G, clusters)
    results['cluster_UB'] = cluster_UB
    
    return results