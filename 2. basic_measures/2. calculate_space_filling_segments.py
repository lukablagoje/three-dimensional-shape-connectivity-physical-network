import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ast
import seaborn as sns
from scipy.stats import linregress
import math
from physical_properties_functions import *
import time


# Initialize a list with the name of the networks to process.
name_list = ['root_1', 'root_2', 'human_neuron', 'zebrafish_neuron', 'monkey_neuron',
             'rat_neuron', 'anthill', 'vascular_1', 'vascular_2', 'vascular_3',
             'mitochondrial', 'fruit_fly_1', 'fruit_fly_2', 'fruit_fly_3', 'fruit_fly_4']

# Define the base path for data files.
path = '../1. data/3. final_data/'

# Initialize a dictionary to hold the final results.
final_results = {}
multiplication_factor = 1
starting_proportion = 5
# Process each network name in the list.
for name in name_list:
    print('**** Network:', name)
    print('Starting computation...')
    final_results[name] = {}

    # Load the paths data for the current network and remove any duplicate rows.
    skeleton_paths = pd.read_csv(path + name + '.paths.csv', index_col=[0])
    skeleton_paths.drop_duplicates(inplace=True)

    # Calculate the bounding box and its volume based on the merged skeleton paths.
    bounds = bounding_box_from_merged_skeleton(skeleton_paths)
    bounding_box_volume = np.prod([abs(bounds[i] - bounds[i+1]) for i in range(0, len(bounds), 2)])
    longest_side = max(abs(bounds[i] - bounds[i+1]) for i in range(0, len(bounds), 2))
    shortest_side = min(abs(bounds[i] - bounds[i+1]) for i in range(0, len(bounds), 2))

    # Calculate distance and volume for each path in the skeleton.
    skeleton_paths['distance'] = np.linalg.norm(skeleton_paths[['x_2', 'y_2', 'z_2']].values - skeleton_paths[['x_1', 'y_1', 'z_1']].values, axis=1)
    skeleton_paths['volume'] = (1/3) * np.pi * (skeleton_paths['radius_1']**2 + skeleton_paths['radius_2']**2 + skeleton_paths['radius_1']*skeleton_paths['radius_2']) * skeleton_paths['distance']

    # Determine the resolution for segments based on their size.
    segment_resolution = min(skeleton_paths['distance'].mean(), skeleton_paths['radius_1'].mean(), skeleton_paths['radius_2'].mean()) / 2
    cube_side = np.max(skeleton_paths['distance'])
    final_results[name]['cube_side'] = cube_side
    final_results[name]['starting_proportion'] = starting_proportion
    final_results[name]['multiplication_factor'] = multiplication_factor
    final_results[name]['segment_resolution_distance'] = segment_resolution

    #multiplied_skeleton['geometry'] = multiplied_skeleton.apply(lambda row: LineString([(row['x_1'], row['y_1'], row['z_1']), (row['x_2'], row['y_2'], row['z_2'])]), axis=1)
    
    # Initialize lists to store computed densities.
    space_filling_list = []
    link_count_list = []
    local_space_filling_greater_count = 0
    
    # Generate ranges for the 3D grid based on bounding box limits and cube side length.
    # Iterate over each cube side value to calculate densities within cubes of that edge length.
    num_cubes = 10
    grid_x = np.linspace(bounds[0], bounds[1], num=num_cubes+1)
    grid_y = np.linspace(bounds[2], bounds[3], num=num_cubes+1)
    grid_z = np.linspace(bounds[4], bounds[5], num=num_cubes+1)
    
    # Compute the 3D grid dimensions
    x_len, y_len, z_len = len(grid_x) - 1, len(grid_y) - 1, len(grid_z) - 1
    

    # Convert grid ranges to a more efficient structure
    grid_ranges = pd.MultiIndex.from_product([range(x_len), range(y_len), range(z_len)], names=['x', 'y', 'z']).to_frame(index=False)

    # Example: Initialize results
    link_count_list = []
    space_filling_list = []
    segment_count_list =[]
    
    # Initialize timing
    start_time = time.time()

    # Loop through each grid cell using vectorized operations where possible
    for index, grid in grid_ranges.iterrows():
        if index%1000 == 0:
            elapsed_time = time.time() - start_time
            print(f'Network {name} progress: {np.round(index / len(grid_ranges), 3)} |',
                  f'Count local_space_filling > 1: {np.round(local_space_filling_greater_count)} |',
                  f'Elapsed time: {np.round(elapsed_time, 2)} seconds')
    
        i, j, k = grid['x'], grid['y'], grid['z']
        # Define the cube boundary
        x_min, x_max =  grid_x[i],  grid_x [i+1]
        y_min, y_max =  grid_y[j],  grid_y [j+1]
        z_min, z_max =  grid_z[k], grid_z[k+1]
        small_box_volume = abs(x_min - x_max)*abs(y_min - y_max) * abs(z_min-z_max)
        # Query to filter segments inside the current cube
        query_filter = f'((x_1 >= {x_min}) & (x_1 < {x_max}) & (x_2 >= {x_min}) & (x_2 < {x_max})) & \
                  ((y_1 >= {y_min}) & (y_1 < {y_max}) & (y_2 >= {y_min}) & (y_2 < {y_max})) & \
                  ((z_1 >= {z_min}) & (z_1 < {z_max}) & (z_2 >= {z_min}) & (z_2 < {z_max}))'

        # Now filter these segments further
        segments_inside =  skeleton_paths.query(query_filter)
        segment_count_list.append(len(segments_inside))

        link_count = 0 if len(segments_inside) == 0 else len(set(segments_inside['path_id'].values))
        link_count_list.append(link_count)

        # Computing link counts and space filling density
        link_count = len(set(segments_inside['path_id'])) if not segments_inside.empty else 0
        local_space_filling = segments_inside['volume'].sum() /  small_box_volume  if not segments_inside.empty else 0
        if local_space_filling > 1 :
            local_space_filling_greater_count += 1
            space_filling_list.append(1)
        else:
            space_filling_list.append(local_space_filling) # Ensuring density doesn't exceed 1
            
    # Store computed density values in the results dictionary.
    final_results[name]['space_filling_list'] = space_filling_list
    final_results[name]['link_count_list'] = link_count_list
    final_results[name]['segment_count_list'] = link_count_list
    final_results[name]['count_local_space_filling_greater'] = local_space_filling_greater_count
    final_results[name]['number_of_cubes'] = len(grid_ranges)
    # Save the density results to a pickle file.
    with open('1. space_filling_results/'+name + "_space_filling_results_segments_n_cubes_" + str( num_cubes)+".pkl", "wb") as h:
        pickle.dump(final_results[name], h)