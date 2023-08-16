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
from networkx import grid_2d_graph
from networkx import grid_graph
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
name_list = ['human_neuron','rat_neuron',
            'monkey_neuron','zebrafish_neuron',
    'vascular_2','vascular_3','vascular_1','mitochondrial','root_1','root_2','anthill','fruit_fly_2','fruit_fly_3','fruit_fly_1','fruit_fly_4']


df_dict = {'name':[],
    'N_number_of_nodes_2D_min': [],
    'N_number_of_nodes_3D_min': [],
    'N_number_of_nodes_2D_max': [],
    'N_number_of_nodes_3D_max': [],
    'N_number_of_nodes': [],
    'log_N': [],
    'network_diameter_1D':[],
    'network_diameter_2D_mean': [],
    'network_diameter_2D_std': [],
    'network_diameter_3D_mean': [],
    'network_diameter_3D_std': []}
path_save = '2. link_degree_betw_results/'
link_degree_name_dict = {}
link_betw_pathid_name_dict = {}
tree_dict = {}
cycle_dict = {}
degree_dict = {}
        
link_network_results = {}
for name in name_list:
    print('******',name)
    df_dict['name'].append(name)
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


    null_model_2D_diameter_list = []
    null_model_3D_diameter_list = []
    trials = 3
    for i in range(0,trials):
        print(name,i,'/',trials)
        if i == 0:
            l_2d = int(np.floor((np.sqrt(len(G.nodes)))))
            l_3d = int(np.floor((len(G.nodes)**(1/3))))
            G_2D = grid_2d_graph(l_2d ,l_2d )
            G_3D = grid_graph(dim=(l_3d, l_3d, l_3d))
            df_dict['N_number_of_nodes_2D_min'].append(len(G_2D.nodes()))
            df_dict['N_number_of_nodes_3D_min'].append(len(G_3D.nodes()))   
        elif i == 1:
            l_2d = int(np.ceil((np.sqrt(len(G.nodes)))))
            l_3d = int(np.ceil((len(G.nodes)**(1/3))))
            G_2D = grid_2d_graph(l_2d ,l_2d )
            G_3D = grid_graph(dim=(l_3d, l_3d, l_3d))      
            df_dict['N_number_of_nodes_2D_max'].append(len(G_2D.nodes()))
            df_dict['N_number_of_nodes_3D_max'].append(len(G_3D.nodes()))
        else:
            G_1D = nx.path_graph(len(G.nodes))
        null_model_2D_diameter = nx.diameter(G_2D)
        null_model_3D_diameter = nx.diameter(G_3D)
        null_model_2D_diameter_list.append(null_model_2D_diameter)
        null_model_3D_diameter_list.append(null_model_3D_diameter)
    null_model_1D_diameter = nx.diameter(G_1D)
    average_null_2D_model_diameter = np.mean(null_model_2D_diameter_list)
    std_null_2D_model_diameter =  np.std(null_model_2D_diameter_list)
    average_null_3D_model_diameter = np.mean(null_model_3D_diameter_list)
    std_null_3D_model_diameter =  np.std(null_model_3D_diameter_list)
    df_dict['N_number_of_nodes'].append(len(G.nodes()))
    df_dict['log_N'].append(np.log(len(G.nodes())))
    df_dict['network_diameter_1D'].append( null_model_1D_diameter)
    df_dict['network_diameter_2D_mean'].append(average_null_2D_model_diameter)
    df_dict['network_diameter_2D_std'].append(std_null_2D_model_diameter)
    df_dict['network_diameter_3D_mean'].append(average_null_3D_model_diameter)
    df_dict['network_diameter_3D_std'].append(std_null_3D_model_diameter)
    df_lattice = pd.DataFrame(df_dict)
    df_lattice.to_csv('network_measures_lattice.csv')

df_lattice = pd.DataFrame(df_dict)
df_lattice.to_csv('network_measures_lattice.csv')