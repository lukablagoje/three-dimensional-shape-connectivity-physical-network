import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
import random
import timeit
from scipy.spatial import distance as dst
import itertools
import pickle
import scipy.stats as stats

plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rc('xtick', labelsize=14) 
plt.rc('ytick', labelsize=14) 
import seaborn as sns
#sns.set_theme(style="whitegrid")
from scipy.stats.stats import spearmanr
    
import networkx as nx
#name_list = ['human_neuron','rat_neuron',
#            'monkey_neuron','zebrafish_neuron',
#    'vascular_2','vascular_3','vascular_1','mitochondrial','root_1','root_2','anthill','fruit_fly_1','fruit_fly_2','fruit_fly_3'] 
#name_list = ['human_neuron','rat_neuron',
#            'monkey_neuron','zebrafish_neuron', 'vascular_2','vascular_3','vascular_1','mitochondrial','root_1','root_2','anthill']
name_list=['fruit_fly_1']
color_dict = {'root_1':'green','root_2':'olive','mitochondrial':'y','vascular_1':'red','vascular_2':'crimson','vascular_3':'salmon',
              'fruit_fly_2':'purple','fruit_fly_3':'pink','tree':'brown','rat_neuron':'magenta','human_neuron':'indigo',
              'anthill':'silver','fruit_fly_1':'deeppink','fruit_fly_4':'plum','zebrafish_neuron':'blue','monkey_neuron':'teal'}
t_std = {}
t_1_distribution_dict = {}

#path_save = '4. directed_metagraph_processed/'
name_corr_df = {}
name_corr_p_df = {}
t_means = {}
t_1_intersections = {}
t_1_intersections_links = {}
t_dict_std = {}
intersections_all_results = {}
for name in name_list:
    print(name)
    path_source = "3. directed_metagraph_results/"
    #infile = open(path_source +  name + "_directed_metagraph_dict_results.pkl",'rb')
    infile = open(name + "_directed_metagraph_dict_results.pkl",'rb')
    intersection_metagraph = pickle.load(infile)
    
    #Creating metagraph
#    directed_metagraph = {}
#    t = 1
#    trial = 0
#    for path_id in intersection_metagraph[t][trial].keys():
#        directed_metagraph[path_id] = {}
#        for path_id_2 in intersection_metagraph[t][trial].keys():
#            directed_metagraph[path_id][path_id_2] = 0
#    path_id_list_dictionary = {}
#    for path_id in intersection_metagraph[t][trial].keys():
#        path_id_list_dictionary[path_id] = []

#    nr_trials = len(intersection_metagraph[t].keys())
#    for trial in range(0,nr_trials):
#        for path_id_shuffled in intersection_metagraph[t][trial].keys():
#            for path_id_hit in intersection_metagraph[t][trial][path_id_shuffled]:
#                 path_id_list_dictionary[path_id_shuffled].append(path_id_hit)
    #Weighted metagraph

#    for path_id_1 in path_id_list_dictionary.keys():
#        for path_id_2 in path_id_list_dictionary[path_id_1]:
#                directed_metagraph[path_id_1][path_id_2] +=1
#    t = 1
#    trial = 0
#    for path_id_1 in intersection_metagraph[t][trial].keys():
#       for path_id_2 in intersection_metagraph[t][trial].keys():
#                directed_metagraph[path_id_1][path_id_2] /= nr_trials
            
    #Out degree measure
#    out_degree_weight_sum = {}
#    for path_id in directed_metagraph.keys():
#       out_degree_weight_sum[path_id] = np.sum(list(directed_metagraph[path_id].values()))
    

    #In degree measure
#    in_degree_weight_sum = {}
#    for path_id_1 in directed_metagraph.keys():
#        for path_id_2 in directed_metagraph.keys():
#                in_degree_weight_sum[path_id_2] = 0 

        
#    for path_id_1 in directed_metagraph.keys():
#        for path_id_2 in directed_metagraph.keys():
#                in_degree_weight_sum[path_id_2] += directed_metagraph[path_id_1][path_id_2] 
    print('Creating metagraph')
    t=1
    in_degree_weight_sum = {}
    out_degree_weight_sum = {}
    for trial in list(intersection_metagraph[t].keys()):
        for path_id_1 in list(intersection_metagraph[t][trial].keys()):
            for path_id_2 in list(intersection_metagraph[t][trial][path_id_1]):
                out_degree_weight_sum[path_id_1] = 0 
                in_degree_weight_sum[path_id_2] = 0 

    for trial in list(intersection_metagraph[t].keys()):
        for path_id_1 in list(intersection_metagraph[t][trial].keys()):
            for path_id_2 in list(intersection_metagraph[t][trial][path_id_1]):
                out_degree_weight_sum[path_id_1] += 1
                in_degree_weight_sum[path_id_2] += 1 


    nr_trials = len(intersection_metagraph[t].keys())
    for path_id in out_degree_weight_sum.keys():
        out_degree_weight_sum[path_id] /= nr_trials

    for path_id in in_degree_weight_sum.keys():
        in_degree_weight_sum[path_id] /= nr_trials
        
    name_network_measures_dict = {}
    print('****** Network Measures',name)
    network_measures_dict = {}
    path_source_2 = '../1. data/3. final_data/'
    path_save = '4. directed_metagraph_weights/'
    # Connectome properties
    link_paths = pd.read_csv(path_source_2 + name+'.paths.csv',index_col=[0])
    path_bodyid_list = link_paths[['path_id','bodyId_pre','bodyId_post']].values.tolist()
    all_body_id_list = list(set(list(set(link_paths['bodyId_pre'].values.tolist())) + list(set(link_paths['bodyId_post'].values.tolist())) ))
    path_bodyid_dict = {}
    for bodyid in all_body_id_list:
        path_bodyid_dict[bodyid] = []
    for path_body in path_bodyid_list:
        path_bodyid_dict[path_body[1]].append(path_body[0])
        path_bodyid_dict[path_body[2]].append(path_body[0])
    for key in path_bodyid_dict.keys():
        path_bodyid_dict[key] = list(set(path_bodyid_dict[key]))

    adjacent_paths = {}
    for path_body in path_bodyid_list:
        path_id_1 = path_body[0]
        path_list_1 = path_bodyid_dict[path_body[1]] 
        path_list_2 = path_bodyid_dict[path_body[2]]
        adjacent_paths[path_body[0]] = {}
        path_id_list = list(set(path_list_1 + path_list_2))
        path_id_list.remove(path_id_1)
        for path_id_2 in path_id_list:
            adjacent_paths[path_id_1][path_id_2] = 1

    #Out degree measure
    connectome_out_degree_weight_sum = {}
    for path_id in adjacent_paths.keys():
        connectome_out_degree_weight_sum[path_id] = np.sum(list(adjacent_paths[path_id].values()))

    connectome_in_degree_weight_sum = {}
    for path_id_1 in adjacent_paths.keys():
        for path_id_2 in adjacent_paths.keys():
            if path_id_1 != path_id_2:
                connectome_in_degree_weight_sum[path_id_2] = 0 

    for path_id_1 in  adjacent_paths.keys():
        for path_id_2 in  adjacent_paths[path_id_1].keys():
            if path_id_1 != path_id_2:
                connectome_in_degree_weight_sum [path_id_2] += adjacent_paths[path_id_1][path_id_2] 
    #average_in_degree = np.mean(list(connectome_in_degree_weight_sum.values()))
    #average_out_degree = np.mean(list(connectome_out_degree_weight_sum.values()))

    #diff_in_degree_weight_sum = {}
    #diff_out_degree_weight_sum = {}
    #for path_id in connectome_out_degree_weight_sum.keys():
     #   diff_out_degree_weight_sum[path_id] = (out_degree_weight_sum[path_id] -  connectome_out_degree_weight_sum[path_id])
     #   diff_in_degree_weight_sum[path_id] = (in_degree_weight_sum[path_id] -  connectome_in_degree_weight_sum[path_id])
    
    weighted_results = {}
    weighted_results['metagraph_out_degree_weight_sum'] = out_degree_weight_sum 
    weighted_results['metagraph_in_degree_weight_sum'] = in_degree_weight_sum 
    weighted_results['connectome_out_degree_weight_sum'] = connectome_out_degree_weight_sum 
    weighted_results['connectome_in_degree_weight_sum'] = connectome_in_degree_weight_sum 
    with open(path_save+ name + '_weighted_results.pickle', 'wb') as handle:
        pickle.dump(weighted_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
