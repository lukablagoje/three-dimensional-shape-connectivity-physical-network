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
name_list = ['human_neuron','rat_neuron',
            'monkey_neuron','zebrafish_neuron',
    'vascular_2','vascular_3','vascular_1','mitochondrial','root_1','root_2','anthill','fruit_fly_2','fruit_fly_3','fruit_fly_1','fruit_fly_4']


df_dict = {'name':[],
           'degree_distribution_mean': [],
 'degree_distribution_std': [],
 '2nd_order_degree_mean': [],
 '2nd_order_degree_std': [],
 'network_diameter': [],
  'average_shortest_path_length':[],
 'average_clustering': [],
 'transitivity': [],
   'N_number_of_nodes':[],
    'log_N':[] ,
    'network_diameter/log(N)':[],
          'max_degree':[],
           'network_diameter_tree_std':[],
           'network_diameter_tree_mean':[],
           'is_not_tree':[]
          }

path_save = '2. link_degree_betw_results/'
link_degree_name_dict = {}
link_betw_pathid_name_dict = {}
tree_dict = {}
cycle_dict = {}
degree_dict = {}
        
link_network_results = {}
for name in name_list:
    print('******',name)
    network_measures_dict = {}
    path_source = '1. network_measures_results/'
    infile = open(path_source+ name + "_network_measures_dict.pkl",'rb')
    network_measures_dict =  pickle.load(infile)
    network_measures_dict['name'] = name
    for key in network_measures_dict.keys():
        if key != 'link_degree_dict' and key!= 'betw_dict' and key !=  'link_degree_dist':
            df_dict[key].append(network_measures_dict[key])
    
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
        G.add_edge(bodyid_edge[0],bodyid_edge[1])
    tree_dict[name] = nx.is_tree(G)
    #print(nx.is_arborescence(G))
    #print(nx.is_branching(G))

    link_degree_dict = {}
    degrees = dict(G.degree())
    degree_dict[name] = degrees
    for path_bodyid in path_bodyid_list:
        link = [path_bodyid[1],path_bodyid[2]]
        pathid = path_bodyid[0]
        link_degree_dict[pathid] = degrees[link[0]] + degrees[link[1]]
                                                              
    link_betw_pathid = {}                                                      
    for pair in network_measures_dict['betw_dict'].keys():
        path_id = bodyid_path_dict[pair]
        link_betw_pathid[path_id] =  network_measures_dict['betw_dict'][pair]
    tree_diameter_list = []
    trials = 20
    for i in range(0,trials):
        print(name,i,'/',trials)
        G_tree =nx.random_tree(len(G.nodes()),seed=i)
        tree_diameter = nx.diameter(G_tree)
        tree_diameter_list.append(tree_diameter)
    average_tree_diameter = np.mean(tree_diameter_list)
    std_tree_diameter =  np.std(tree_diameter_list)
    link_degree_name_dict[name] = link_degree_dict
    link_betw_pathid_name_dict[name] = link_betw_pathid
    #cycle_dict[name] 
    #if nx.is_tree(G) == False:
    #cycle_dict[name]= nx.cycle_basis(G)
    #df_dict['cycle_length_distribution'].append(cycle_lengths)
    df_dict['N_number_of_nodes'].append(len(G.nodes()))
    df_dict['log_N'].append(np.log(len(G.nodes())))
    df_dict['network_diameter/log(N)'].append(network_measures_dict['network_diameter']/np.log(len(G.nodes())))
    df_dict['network_diameter_tree_mean'].append(average_tree_diameter)
    df_dict['network_diameter_tree_std'].append(std_tree_diameter)
   # df_dict['standard_link_ratio'].append( standard_link_ratio)
    #df_dict['non_standard_link_ratio'].append( non_standard_link_ratio )
    if tree_dict[name]:
        df_dict['is_not_tree'].append(0)
    else:
        df_dict['is_not_tree'].append(1)
        
link_network_results['link_degree']  =link_degree_name_dict 
link_network_results['link_betw'] = link_betw_pathid_name_dict
with open("1. network_measures_results/link_degree_betw.pkl", "wb") as h:
     pickle.dump(link_network_results, h)

df = pd.DataFrame(df_dict).sort_values(by=['average_shortest_path_length'])
df['z-score_n_dia'] = (df['network_diameter']-df['network_diameter_tree_mean'])/df_dict['network_diameter_tree_std']

df.to_csv('network_measures.csv')