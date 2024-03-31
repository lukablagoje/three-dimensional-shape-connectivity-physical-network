import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ast
import seaborn as sns
from scipy.stats import linregress
import math
from scipy.spatial import distance as dst
from scipy.spatial import ConvexHull

from physical_properties_functions import *

# Define a list of network names for processing.
name_list = ['root_1', 'root_2', 'human_neuron', 'zebrafish_neuron', 'monkey_neuron',
             'rat_neuron', 'anthill', 'vascular_1', 'vascular_2', 'vascular_3',
             'mitochondrial', 'fruit_fly_1', 'fruit_fly_2', 'fruit_fly_3', 'fruit_fly_4']

# Initialize a dictionary to store final results.
final_results = {}

# Loop through each network name for processing.
for name in name_list:
    print('**** Network:', name)
    print('Starting computation')
    final_results[name] = {}  # Initialize a sub-dictionary for each network.
    path = '../1. data/3. final_data/'  # Define the base path for data files.
    
    # Load skeleton paths from CSV, setting the first column as the index.
    skeleton_paths = pd.read_csv(path + name + '.paths.csv', index_col=[0])
    skeleton_paths.drop_duplicates(inplace=True)  # Remove duplicate rows.

    # Calculate bounding box, volume, and side lengths for the skeleton.
    bounds = bounding_box_from_merged_skeleton(skeleton_paths)
    bounding_box_volume = np.prod([abs(bounds[i] - bounds[i+1]) for i in range(0, len(bounds), 2)])
    longest_side = max(abs(bounds[i] - bounds[i+1]) for i in range(0, len(bounds), 2))
    shortest_side = min(abs(bounds[i] - bounds[i+1]) for i in range(0, len(bounds), 2))

    # Calculate the distance between points in each path.
    skeleton_paths['distance'] = np.sqrt((skeleton_paths['x_2'] - skeleton_paths['x_1'])**2 +  
                                         (skeleton_paths['y_2'] - skeleton_paths['y_1'])**2 +
                                         (skeleton_paths['z_2'] - skeleton_paths['z_1'])**2)

    # Define cube side lengths as a proportion of the shortest side.
    cube_side = shortest_side * 0.32
    cube_side_list = [proportion / 100 * cube_side for proportion in np.arange(15, 85, 5)]

    # Determine the resolution for segment lengths.
    segment_resolution = np.min([skeleton_paths['distance'].mean(),
                                 skeleton_paths['radius_1'].mean(),
                                 skeleton_paths['radius_2'].mean()]) / 2
    final_results[name]['cube_side_list'] = cube_side_list  # Store the list of cube sides.

    print('Skeleton multiplication')
    multiplied_skeleton_part_merged = []  # Initialize list for multiplied skeleton parts.
    
    # Loop through each row in the skeleton paths dataframe.
    for row in skeleton_paths[['pt_id_1', 'x_1', 'y_1', 'z_1', 'radius_1', 'pt_id_2', 'x_2', 'y_2', 'z_2', 'radius_2']].values:
        # Extract start and end points, and their radii.
        start_point = row[1:4]
        radius_start = row[4]
        end_point = row[6:9]
        radius_end = row[9]

        # Calculate distance and ratio for determining multiplication factor.
        point_to_point_dist = dst.euclidean(start_point, end_point)
        ratio = (point_to_point_dist + radius_start + radius_end) / cube_side
        
        # Create a new row with start and end points, including their radii.
        new_row = [list(start_point) + [radius_start] + list(end_point) + [radius_end]]
        multiplied_skeleton_part_merged += new_row  # Add the new row to the list.

    # Convert the list of multiplied parts into a dataframe.
    multiplied_skeleton = pd.DataFrame(multiplied_skeleton_part_merged,
                                       columns=['x_1', 'y_1', 'z_1', 'radius_1', 'x_2', 'y_2', 'z_2', 'radius_2'])
    skeleton_paths = multiplied_skeleton  # Replace the original skeleton paths with the multiplied version.
    print('Skeleton multiplication complete')
    
    # Initialize a list to store point cloud data.
    point_cloud_list = []
    # Loop through each row in the modified skeleton paths to create point clouds.
    for i, row in enumerate(skeleton_paths[['x_1', 'y_1', 'z_1', 'radius_1', 'x_2', 'y_2', 'z_2', 'radius_2']].values):
        print('**** Network:', name)
        print('Computing point clouds', i, len(skeleton_paths))
        # Generate point clouds for each segment in the skeleton.
        point_cloud_list += create_thick_line_point_list(row[0:3], row[4:7], segment_resolution, row[3], row[7])
    
    # Convert the list of point clouds into a dataframe.
    point_cloud_skeleton_paths = pd.DataFrame(point_cloud_list, columns=['x', 'y', 'z'])
    
    filled_boxes_count_list = []  # Initialize list to store counts of filled boxes.
    # Loop through each cube side to calculate filled box counts.
    for cube_side in cube_side_list:
        print('**** Network:', name, 'Cube side', cube_side, 'out of', cube_side_list)
        print('Cube side ratio:', cube_side / longest_side)

        # Generate ranges for the 3D grid based on bounding box limits and cube side length.
        x_axis = np.arange(bounds[0], bounds[1] + cube_side, cube_side)
        y_axis = np.arange(bounds[2], bounds[3] + cube_side, cube_side)
        z_axis = np.arange(bounds[4], bounds[5] + cube_side, cube_side)

        filled_boxes_count = 0  # Initialize count of filled boxes.
        empty_boxes_count = 0  # Initialize count of empty boxes.
        
        # Loop through each cube in the 3D grid to count filled and empty boxes.
        for i in range(len(x_axis) - 1):
            for j in range(len(y_axis) - 1):
                for k in range(len(z_axis) - 1):
                    # Identify points within the current cube.
                    points_inside = point_cloud_skeleton_paths[
                        (point_cloud_skeleton_paths['x'] >= x_axis[i]) &
                        (point_cloud_skeleton_paths['x'] <= x_axis[i + 1]) &
                        (point_cloud_skeleton_paths['y'] >= y_axis[j]) &
                        (point_cloud_skeleton_paths['y'] <= y_axis[j + 1]) &
                        (point_cloud_skeleton_paths['z'] >= z_axis[k]) &
                        (point_cloud_skeleton_paths['z'] <= z_axis[k + 1])
                    ]

                    # Update counts based on whether the current cube contains points.
                    if len(points_inside) == 0:
                        empty_boxes_count += 1
                    else:
                        filled_boxes_count += 1
        
        filled_boxes_count_list.append(filled_boxes_count)  # Store the count of filled boxes for each cube side.

    # Store the list of filled box counts in the final results.
    final_results[name]['filled_boxes_count_list'] = filled_boxes_count_list

    # Save the computed results to a pickle file.
    with open("2. fractal_dimension_results/" + name + "_basic_results.pkl", "wb") as h:
           pickle.dump(final_results[name], h)    