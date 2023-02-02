import networkx as nx
import geopandas as gpd
from util import get_state
import matplotlib.pyplot as plt

def draw_clusters(G, clusters, sizes, png_filename):

    # G corresponds to which state?
    node = nx.utils.arbitrary_element(G.nodes)
    fips = G.nodes[node]['GEOID20'][0:2]
    state = get_state(fips)
    
    # get county-level geopandas dataframe 
    shp_filepath = 'C:\\districting-data-2020\\'
    shp_filename = state + "_county.shp"
    df = gpd.read_file( shp_filepath + shp_filename )
    
    # labeling[i]=j if node i is assigned to cluster j
    labeling = { i : j for j in range(len(clusters)) for i in clusters[j] }
    node_with_this_geoid = { G.nodes[i]["GEOID20"] : i for i in G.nodes }
    color = [ -1 for u in G.nodes ]
    
    # find a coloring of the shrunk cluster graph
    shrunk_graph = nx.Graph()
    shrunk_graph.add_nodes_from( range(len(sizes) ) )
    shrunk_graph.add_edges_from( list( { (labeling[i],labeling[j]) for i,j in G.edges if labeling[i]!=labeling[j] } ) )
    
    # we actually use a coloring of the square of the shrunk graph.
    # this avoids issues with the same color given to different clusters that meet at a point.
    square_graph = nx.Graph()
    square_graph.add_nodes_from( range(len(sizes) ) )
    for i in square_graph.nodes:
        dist = nx.shortest_path_length(shrunk_graph, source=i)
        square_graph.add_edges_from( [ (i,j) for j in dist.keys() if j > i and dist[j] <= 2 ] )
    
    cluster_color = nx.greedy_color(square_graph, strategy='smallest_last')

    # pick a position u in the dataframe
    for u in range(len(G.nodes)):
        geoid = df['GEOID20'][u]

        # what node in G has the same geoid?
        i = node_with_this_geoid[geoid]

        # position u in the dataframe should be given 
        #  the same district label as i has in labeling
        color[u] = cluster_color[ labeling[i] ]

    # Now add the assignments to a column of the dataframe and map it
    df['color'] = color
    my_fig = df.plot(column='color',edgecolor='white',cmap='Paired').get_figure()
    RESIZE_FACTOR = 3
    my_fig.set_size_inches(my_fig.get_size_inches()*RESIZE_FACTOR)
    plt.axis('off')
    plt.tight_layout()

    # show the size of each cluster on the map
    for p in range(len(clusters)):
        cluster = clusters[p]
        #denom = len(cluster)
        denom = sum( G.nodes[i]['area'] for i in cluster )
        x = sum( G.nodes[i]['C_X']*G.nodes[i]['area'] for i in cluster ) / denom
        y = sum( G.nodes[i]['C_Y']*G.nodes[i]['area'] for i in cluster ) / denom
        plt.annotate( str(sizes[p]), xy=(x,y), size=20, horizontalalignment='center', verticalalignment='center' )
    
    my_fig.savefig(png_filename)
    return
    