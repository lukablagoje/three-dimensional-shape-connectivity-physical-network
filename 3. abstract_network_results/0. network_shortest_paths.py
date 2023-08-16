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
name_list = ['human_neuron','rat_neuron',
            'monkey_neuron','zebrafish_neuron',
    'vascular_2','vascular_3','vascular_1','mitochondrial','anthill','root_1','root_2','fruit_fly_2','fruit_fly_3','fruit_fly_1','fruit_fly_4'] 
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
    link_degree = {}
    for link in G.edges():
        link_degree[link] = degrees[link[0]] + degrees[link[1]]
    network_measures_dict['shortest_paths'] = betw_dict    
    with open("2. shortest_paths_results/"+ name +"_shortest_paths_dict.pkl", "wb") as h:
        pickle.dump(network_measures_dict, h)