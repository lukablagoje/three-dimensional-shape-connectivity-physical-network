import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
import math
import os
import timeit
import random
from physical_properties_functions import *


def randomize_segments_with_thickness(segment_list, radius_list, seed):
    """
    Randomizes the order of segments while preserving their thickness,
    simulating variations in segment arrangements without altering their physical properties.
    
    Parameters:
    - segment_list (list): A list of segments where each segment is represented as a pair of points (start and end).
    - radius_list (list): A list of radii corresponding to each segment in segment_list.
    - seed (int): A seed for the random number generator to ensure reproducibility.
    
    Returns:
    - randomized_segment_list (list): A new list of segments that have been re-ordered based on the randomization process.
    - radius_list (tuple): A tuple of radii reordered in accordance with the randomized segments.
    """
    # Initialize list to store directional vectors of segments
    diff_vector_list = [np.array(segment[1]) - np.array(segment[0]) for segment in segment_list]
    
    # Set the random seed for reproducibility
    np.random.seed(seed)
    
    # Pair each directional vector with its corresponding radius
    c = list(zip(diff_vector_list, radius_list))
    
    # Randomly shuffle the pairs to randomize the order of segments
    np.random.shuffle(c)
    
    # Unzip the shuffled pairs back into directional vectors and radii
    diff_vector_list, radius_list = zip(*c)
    
    # Reconstruct the randomized segments from the shuffled directional vectors
    randomized_segment_list = []
    starting_point = segment_list[0][0]
    for diff_vector in diff_vector_list:
        new_second_point = starting_point + diff_vector
        randomized_segment_list.append([np.array(starting_point), np.array(new_second_point)])
        starting_point = new_second_point
        
    return randomized_segment_list, radius_list

def intersecting_bodyid_from_kd_tree(intersections, dict_of_all_points):
    """
    Identifies unique bodies (or segments) intersecting based on indices from KD-tree queries.
    
    Parameters:
    - intersections (list): A list of lists where each sublist contains indices of points intersecting a query.
    - dict_of_all_points (dict): A dictionary mapping the index of a point in a KD-tree to its corresponding body ID.
    
    Returns:
    - bodyid_intersected (list): A list of unique body IDs that intersect the query.
    """
    # Flatten the list of intersections and remove duplicates to get unique point indices
    point_indices = set(num for sublist in intersections for num in sublist)
    
    # Initialize a list to store the body IDs of intersecting points
    bodyid_intersected = []
    for points_index in point_indices:
        # Retrieve the body ID corresponding to each intersecting point index
        for i, point_key in enumerate(sorted(dict_of_all_points.keys())):
            if point_key >= points_index:
                bodyid_intersected.append(dict_of_all_points[point_key])
                break
                
    return bodyid_intersected

def list_dict_all_points(points_bodyid):
    """
    Creates a list of all points and a dictionary mapping point indices to their respective body IDs.
    
    Parameters:
    - points_bodyid (dict): A dictionary where keys are body IDs and values are lists of points associated with each body.
    
    Returns:
    - list_of_all_points (list): A flattened list containing all points from all bodies.
    - dict_of_all_points (dict): A dictionary mapping the index of each point in list_of_all_points to its corresponding body ID.
    """
    list_of_all_points = []
    dict_of_all_points = {}
    current_index = 0
    
    for key, points in points_bodyid.items():
        list_of_all_points += points
        for _ in points:
            # Map the end index of the current body's points to the body ID
            dict_of_all_points[current_index] = key
            current_index += 1
            
    return list_of_all_points, dict_of_all_points
# Define paths for source data and where to save results
path_save = '1. directed_metagraph_results/'
path_source = '../1. data/3. final_data/'

# List of names for the networks to be analyzed
name_list = ['human_neuron', 'rat_neuron', 'monkey_neuron', 'zebrafish_neuron', 'vascular_2', 'vascular_3', 'vascular_1', 'mitochondrial', 'root_1', 'root_2', 'anthill', 'fruit_fly_2', 'fruit_fly_3', 'fruit_fly_4']

# Number of trials for simulation
N_trials = 5

# Loop through each network to process its data
for name in name_list:
    print(name)
    # Load the skeleton paths data
    skeleton_paths = pd.read_csv(path_source + name + '.paths.csv', index_col=0)
    # Calculate the Euclidean distance between path segments
    skeleton_paths['distance'] = np.sqrt((skeleton_paths['x_2'] - skeleton_paths['x_1'])**2 + (skeleton_paths['y_2'] - skeleton_paths['y_1'])**2 + (skeleton_paths['z_2'] - skeleton_paths['z_1'])**2)
    # Remove duplicate paths
    skeleton_paths.drop_duplicates(inplace=True)

    # Calculate mean radii and segment distances
    radius_mean_1 = skeleton_paths['radius_1'].mean()
    radius_mean_2 = skeleton_paths['radius_2'].mean()
    segment_distance_mean = skeleton_paths['distance'].mean()
    segment_resolution = np.min([segment_distance_mean, radius_mean_1, radius_mean_2])

    # Initialize dictionaries for segments and radii paths
    segments_path_dict = {}
    radius_path_dict = {}
    path_id_list = list(set(skeleton_paths['path_id'].values.tolist()))

    # Populate segment and radius dictionaries with path data
    for path_id in path_id_list:
        segment_list = skeleton_paths[skeleton_paths['path_id'] == path_id][['x_1', 'y_1', 'z_1', 'x_2', 'y_2', 'z_2']].values.tolist()
        radius_list = skeleton_paths[skeleton_paths['path_id'] == path_id][['radius_1', 'radius_2']].values.tolist()
        segments_path_dict[path_id] = [[segment[0:3], segment[3:]] for segment in segment_list]
        radius_path_dict[path_id] = radius_list

    # Create original KD-trees for each path
    kd_tree_original = {}
    for path_id in path_id_list:
        print('Creating original KD-trees', path_id)
        thick_line_list = []
        for segment, radius in zip(segments_path_dict[path_id], radius_path_dict[path_id]):
            if not np.array_equal(segment[0], segment[1]):
                thick_line_list += create_thick_line_point_list(segment[0], segment[1], segment_resolution, radius[0], radius[1])
        kd_tree_original[path_id] = thick_line_list

    # Define adjacency relations between paths
    adjacent_paths = {}
    for path_body in skeleton_paths[['path_id', 'source', 'target']].values.tolist():
        adjacent_paths[path_body[0]] = list(set(path_bodyid_dict[path_body[1]] + path_bodyid_dict[path_body[2]]) - {path_body[0]})


    t_list = [1]  # List of thickness multipliers to test
    start_trial = 0  # Starting index for trials, allows for resuming incomplete simulations

    # Iterate over each thickness multiplier in the test list
    for t in t_list:
        t_dict_results[t] = {}  # Initialize a dictionary to store results for the current thickness multiplier
        # Iterate over each path ID in the network to construct KD-trees representing the spatial distribution of paths
        for count, path_id in enumerate(path_id_list):
            print(count, len(path_id_list))
            segment_list = segments_path_dict[path_id]  # Retrieve segment coordinates for the current path
            radius_list = radius_path_dict[path_id]  # Retrieve radii for the segments in the current path

            # Generate a list of points representing a "thickened" version of the path
            thick_line_list = []
            for i, segment in enumerate(segment_list):
                if not np.array_equal(segment[0], segment[1]):  # Ensure the segment is not a point
                    # Increase the radius by the thickness multiplier and generate points along the segment
                    thick_line_list += create_thick_line_point_list(segment[0], segment[1], segment_resolution, radius_list[i][0] * t, radius_list[i][1] * t)
            kd_tree_original[path_id] = thick_line_list  # Store the generated point list in a KD-tree

        # Combine all points from all paths into a single KD-tree for efficient global spatial queries
        list_of_all_points, dict_of_all_points = list_dict_all_points(kd_tree_original)
        kd_tree_large = KDTree(list_of_all_points, balanced_tree=True, compact_nodes=True)

        # Perform spatial queries for each trial to identify path intersections
        for trial in range(start_trial, N_trials):
            start_time = timeit.default_timer()  # Record the start time of the trial for performance monitoring
            kd_tree_pathid = {}  # Initialize a dictionary to store KD-trees for each path in the current trial

            # Randomize segments and radii for each path to simulate variability
            for count, path_id in enumerate(path_id_list):
                print('Point cloud, trial:', trial, 'path count:', count, len(path_id_list))
                segment_list, radius_list = randomize_segments_with_thickness(segments_path_dict[path_id], radius_path_dict[path_id], trial * (path_id + 1) + count)

                # Generate a thickened point cloud for the randomized path
                thick_line_list = []
                for i, segment in enumerate(segment_list):
                    if not np.array_equal(segment[0], segment[1]):
                        thick_line_list += create_thick_line_point_list(segment[0], segment[1], segment_resolution, t * radius_list[i][0], t * radius_list[i][1])
                kd_tree_pathid[path_id] = KDTree(thick_line_list, balanced_tree=True, compact_nodes=True)

            print('KD preparation for trial', trial)
            kd_tree_radius = np.min([segment_distance_mean, t * radius_mean_1, t * radius_mean_2])

            # Execute spatial queries to detect intersections among paths
            intersections_paths = {}
            for count, path_id in enumerate(path_id_list):
                print('Distance computation, trial:', trial, 'path count:', count, len(path_id_list))
                intersections = kd_tree_pathid[path_id].query_ball_tree(kd_tree_large, r=kd_tree_radius)
                bodyid_intersected = intersecting_bodyid_from_kd_tree(intersections, dict_of_all_points)
                final_bodyid_list = list(set(bodyid_intersected + adjacent_paths[path_id]))
                final_bodyid_list.remove(path_id)  # Exclude the current path from its own intersection list
                intersections_paths[path_id] = final_bodyid_list

            t_dict_results[t][trial] = intersections_paths  # Store the intersections result for the current trial
            print('Parameter t', t, 'time for trial', trial, timeit.default_timer() - start_time)

        # Save the results for all trials at the current thickness multiplier
        with open(path_save + name + '_directed_metagraph_dict_results.pkl', "wb") as h:
            pickle.dump(t_dict_results, h)