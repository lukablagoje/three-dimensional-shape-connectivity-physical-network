import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from networkx.algorithms.connectivity import local_edge_connectivity
import networkx as nx
import seaborn as sns    
import scipy.stats as stats
import pandas as pd
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
#name_list = ['human_neuron','rat_neuron',
#            'monkey_neuron','zebrafish_neuron',
#    'vascular_2','vascular_3','vascular_1','mitochondrial','anthill','root_1','root_2','fruit_fly_2','fruit_fly_3','fruit_fly_1','fruit_fly_4']
#name_list = ['fruit_fly_4']
name_list =['fruit_fly_2','fruit_fly_3','fruit_fly_1','fruit_fly_4']
for name in name_list:
    print('******',name)
    network_measures_dict = {}
    path_source = '../1. data/3. final_data/'
    link_paths = pd.read_csv(path_source +name + '.paths.csv',index_col=[0])
    path_bodyid_list = link_paths[['path_id','bodyId_pre','bodyId_post']].drop_duplicates().values.tolist()
    G = nx.Graph()
    path_bodyid_dict = {}  
    bodyid_path_dict = {}
    for path_bodyid in path_bodyid_list:
        path_bodyid_dict[path_bodyid[0]] =  (path_bodyid[1],path_bodyid[2])
        bodyid_path_dict[(path_bodyid[1],path_bodyid[2])] = path_bodyid[0]
        bodyid_path_dict[(path_bodyid[2],path_bodyid[1])] = path_bodyid[0] 
    for path_id,bodyid_edge in path_bodyid_dict.items():
        #print(path_id,bodyid_edge)
        G.add_edge(bodyid_edge[0],bodyid_edge[1])
    print('Calculating degrees')
    #beetweenness = nx.edge_betweenness_centrality(G,normalized=True)
    degrees =  dict(nx.degree(G))
    #degree_2nd_order_pathid = {}
    print('Calculating 2nd order degrees')
    degrees_2nd_order = {}
    for node in G.nodes():
        degrees_2nd_order[node] =  len(nx.single_source_shortest_path_length(G, node, cutoff=2)) - 1
    link_degree = {}
    for link in G.edges():
        link_degree[link] = degrees[link[0]] + degrees[link[1]]
    #for path_id,bodyid_edge in path_bodyid_dict.items():
    #    neighborhood_nodes = list(set(list(G.neighbors(bodyid_edge[0]))+list(G.neighbors(bodyid_edge[1]))))
    #    total_2nd_order_degree = 0
    #    for node in neighborhood_nodes:
     #       total_2nd_order_degree+= degrees[node]
    print('Calculating network diameter')
    network_diameter = nx.diameter(G)
    print('Calculating average clustering')
    average_clustering = nx.average_clustering(G)
    print('Calculating 2nd order degrees')
    transitivity = nx.transitivity(G)
    print('Average shortest path')
    avg_shortest_path = nx.average_shortest_path_length(G)
    print('Edge betweennees centrality')
    betw_dict = nx.edge_betweenness_centrality(G)
    #network_measures_dict['beetweenneess'] = beetweenness 
    network_measures_dict['degree_distribution_mean'] = np.round(np.mean(list(degrees.values())),2)
    network_measures_dict['degree_distribution_std'] = np.round(np.std(list(degrees.values())),2)
    network_measures_dict['2nd_order_degree_mean'] =   np.round(np.mean(list(degrees_2nd_order.values())),2)
    network_measures_dict['2nd_order_degree_std'] =   np.round(np.std(list(degrees_2nd_order.values())),2)
    network_measures_dict['network_diameter'] =  network_diameter
    network_measures_dict['average_clustering'] =  average_clustering 
    network_measures_dict['transitivity'] =  transitivity
    network_measures_dict['max_degree'] = np.max(list(degrees.values()))
    network_measures_dict['link_degree_dict'] = link_degree
    network_measures_dict['average_shortest_path_length'] = avg_shortest_path
    network_measures_dict['betw_dict'] = betw_dict    
    with open("1. network_measures_results/"+ name +"_network_measures_dict.pkl", "wb") as h:
        pickle.dump(network_measures_dict, h)