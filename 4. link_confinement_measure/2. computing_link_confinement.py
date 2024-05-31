import pandas as pd
import numpy as np
import pickle

# Define a list of names for networks to analyze.
name_list = ['fruit_fly_1']

# Specify the path where the results will be saved.
path_source = "1. collision_results/"
path_save = '2. link_confinement_measure/'

# Loop through each network name to process.
for name in name_list:
    print(name)
    
    # Load the directed metagraph results for the current network.
    infile = open(path_source + name + "_directed_metagraph_dict_results.pkl", 'rb')
    intersection_metagraph = pickle.load(infile)

    print('Creating metagraph')
    
    # Initialize dictionaries to calculate in-degree and out-degree weights.
    in_degree_weight_sum = {}
    out_degree_weight_sum = {}
    
    # Calculate out-degree and in-degree weights for each path in the metagraph.

    for trial, paths in  intersection_metagraph.items():
        for path_id_1, connected_paths in paths.items():
            out_degree_weight_sum[path_id_1] = out_degree_weight_sum.get(path_id_1, 0) + len(connected_paths)
            for path_id_2 in connected_paths:
                in_degree_weight_sum[path_id_2] = in_degree_weight_sum.get(path_id_2, 0) + 1

    # Normalize the in-degree and out-degree weights by the number of trials.
    nr_trials = len(intersection_metagraph.keys())
    for path_id in out_degree_weight_sum:
        out_degree_weight_sum[path_id] /= nr_trials
    for path_id in in_degree_weight_sum:
        in_degree_weight_sum[path_id] /= nr_trials
    
    print('****** Network Measures', name)
    
    # Load connectome properties from CSV file.
    path_source_2 = '../1. data/3. final_data/'
    link_paths = pd.read_csv(path_source_2 + name + '.paths.csv', index_col=[0])
    path_bodyid_list = link_paths[['path_id', 'source', 'target']].values.tolist()

    # Create dictionaries to store adjacent paths for each path id.
    adjacent_paths = {}
    for path_id, source, target in path_bodyid_list:
        adjacent_paths.setdefault(path_id, {})
        for adj_path_id, adj_source, adj_target in path_bodyid_list:
            if path_id != adj_path_id and (source == adj_source or source == adj_target or target == adj_source or target == adj_target):
                adjacent_paths[path_id][adj_path_id] = 1
    
    # Calculate out-degree and in-degree weights for the connectome.
    connectome_out_degree_weight_sum = {path_id: sum(adjacent_paths[path_id].values()) for path_id in adjacent_paths}
    connectome_in_degree_weight_sum = {path_id: 0 for path_id in adjacent_paths}
    for path_id_1 in adjacent_paths:
        for path_id_2 in adjacent_paths[path_id_1]:
            connectome_in_degree_weight_sum[path_id_2] += 1
    
    # Store the weighted results in a dictionary.
    weighted_results = {
        'metagraph_out_degree_weight_sum': out_degree_weight_sum,
        'metagraph_in_degree_weight_sum': in_degree_weight_sum,
        'connectome_out_degree_weight_sum': connectome_out_degree_weight_sum,
        'connectome_in_degree_weight_sum': connectome_in_degree_weight_sum
    }
    
    # Save the weighted results to a pickle file.
    with open(path_save + name + '_weighted_results.pickle', 'wb') as handle:
        pickle.dump(weighted_results, handle, protocol=pickle.HIGHEST_PROTOCOL)