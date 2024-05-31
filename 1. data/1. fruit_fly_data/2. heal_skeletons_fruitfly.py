import numpy as np
import pandas as pd
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import pickle
import timeit
import itertools
from neuprint import fetch_roi_hierarchy
import networkx as nx
from scipy.spatial import distance

# Dictionary mapping each skeleton file to a specific identifier
skeletons_dict = {1: 'POC_96', 4: 'GF(R)_10', 3: 'mALT(L)_20', 2: 'GC_28'}
path_source = '1. neuron_regions_skeleton/'  # Source directory for CSV files
path_save = ''  # Uncomment to specify the save directory

# Iterating through each item in the dictionary
for df_index, skeleton_region in skeletons_dict.items():
    # Load the CSV file for the current region into a DataFrame
    df = pd.read_csv(path_source + skeleton_region + '.csv', index_col=[0])
    
    # Processing segment data from the DataFrame
    segments = df[['link', 'rowId', 'bodyId']].copy()
    segments = segments.rename(columns={'link': 'pt_id1', 'rowId': 'pt_id2'})
    segments = segments[(segments['pt_id1'] != -1) & (segments['pt_id2'] != -1)]
    segments.reset_index(inplace=True)
    segments = segments.rename(columns={'level_0': 'seg_id'})
    segments.drop(['index'], axis=1, inplace=True)
    segments['segment_healing_done'] = 'original'
    # Processing point data from the DataFrame
    points = df[['rowId', 'x', 'y', 'z', 'radius', 'bodyId']].copy()
    points = points.rename(columns={'rowId': 'pt_id'})
    points['diameter'] = points['radius'] * 2
    points.reset_index(inplace=True)

    # Initializing lists for storing processed data
    healed_list_segments = []
    healed_list_points = []

    # Iterating over unique body IDs to process each neuron
    for neuron_index, bodyid in enumerate(set(segments['bodyId'].values)):
        print('Neuron number', neuron_index, 'out of', len(set(segments['bodyId'].values)))
        segments_df = segments[segments['bodyId'] == bodyid].copy()
        points_df = points[points['bodyId'] == bodyid].copy()
        disconnected_skeleton = True
        iteration_step = 0
        
        # Loop to connect components until the skeleton is fully connected
        while disconnected_skeleton:
            iteration_step += 1
            print('Iteration step', iteration_step)
            
            # Create a graph from the segments and check its connectivity
            G = nx.Graph()
            G.add_edges_from(segments_df[['pt_id1', 'pt_id2']].values)
            components = list(nx.connected_components(G))
            print('The number of disconnected components is', len(components))
            if len(components) == 1:
                disconnected_skeleton = False
                break

            # Further processing if there are disconnected components
            node_coordinate_dict = {}
            component_dict = {}
            # Loop through each connected component identified in the graph
            for i, component in enumerate(components):
                # Store the current component in a dictionary with its index as the key
                component_dict[i] = component
                # Iterate over each node in the current component
                for node in component:
                    # Retrieve the x, y, z coordinates for the current node from the dataframe
                    node_coordinate_dict[node] =  points_df.loc[points_df['pt_id'] == node, ['x', 'y', 'z']].values[0]
    
            # Dictionary to track the closest distance between nodes of different components
            component_closest_dict = {}
            # Dictionary to store the actual node pairs with the minimum distance for each component pair
            node_pairs_min_dist = {}
            # Generate all unique pairs of components to compare
            for component_1, component_2 in itertools.combinations(component_dict, 2):
                # Get the coordinates for all nodes in the first component
                nodes_1 = [node_coordinate_dict[node] for node in component_dict[component_1]]
                # Get the coordinates for all nodes in the second component
                nodes_2 = [node_coordinate_dict[node] for node in component_dict[component_2]]
                # Calculate the Euclidean distance between each pair of nodes from the two components
                distances = distance.cdist(nodes_1, nodes_2, 'euclidean')
                # Find the indices of the minimum distance in the distance matrix
                min_index = np.unravel_index(np.argmin(distances), distances.shape)
                # Record the minimum distance for this component pair
                component_closest_dict[(component_1, component_2)] = np.min(distances)
                # Record the specific nodes between which this minimum distance occurs
                node_pairs_min_dist[(component_1, component_2)] = (list(component_dict[component_1])[min_index[0]], list(component_dict[component_2])[min_index[1]])

            # Create segments to connect the closest nodes between components
            if component_closest_dict:  # Ensure there are entries to process
                global_min_distance = np.min(list(component_closest_dict.values()))
                # Find all component pairs with the global minimum distance
                component_pairs_to_heal = [comp_pair for comp_pair, dist in component_closest_dict.items() if dist == global_min_distance]
                # Prepare a new index for adding new rows safely even if the DataFrame is empty
                new_index = max(segments_df.index, default=-1) + 1

                for comp_pair,dist in component_closest_dict.items():
                    if dist == global_min_distance:
                        # Get the node pair that represents the closest points between the two components
                        node_pair = node_pairs_min_dist[comp_pair]
                        # Create a new DataFrame row for the segment to heal
                        new_row = pd.DataFrame({
                            'pt_id1': [node_pair[0]], 
                            'pt_id2': [node_pair[1]], 
                            'segment_healing_done': ['healed']
                        }, index=[new_index])
                        new_index += 1  # Increment the index for the next new row
                        # Append the new row to the segments DataFrame
                        segments_df = pd.concat([segments_df, new_row])
                        break
        print('Healing complete for neuron number', neuron_index, 'out of', len(set(segments['bodyId'].values)))
        segments_df['bodyId'] = bodyid
        segments_df['seg_id'] = np.arange(len(segments_df))
        segments_df.set_index('seg_id', inplace=True)
        segments_df.drop_duplicates(inplace=True)
        segments_df['pt_id1'] = segments_df['pt_id1'].astype(str) + '_' + segments_df['bodyId'].astype(str)
        segments_df['pt_id2'] = segments_df['pt_id2'].astype(str) + '_' + segments_df['bodyId'].astype(str)
        points_df['pt_id'] = points_df['pt_id'].astype(str) + '_' + points_df['bodyId'].astype(str)
        points_df.set_index('pt_id', inplace=True)
        healed_list_segments.append(segments_df)
        healed_list_points.append(points_df)

    # Combine all processed segments and points, then save to CSV
    full_skeleton_segments = pd.concat(healed_list_segments)
    full_skeleton_segments.reset_index(inplace=True)
    full_skeleton_segments.to_csv(path_save + 'fruit_fly_' + str(df_index) + '.segments.csv')
    full_skeleton_points = pd.concat(healed_list_points)
    full_skeleton_points.to_csv(path_save + 'fruit_fly_' + str(df_index) + '.points.csv')