import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
import ast
import seaborn as sns
from scipy.stats import linregress
import math
from scipy.spatial import distance as dst
from scipy.spatial import ConvexHull
from physical_properties_functions import *


# Initialize a list with the name of the networks to process.
name_list = ['fruit_fly_3']

# Define the base path for data files.
path = '../1. data/3. final_data/'

# Initialize a dictionary to hold the final results.
final_results = {}

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
    
    # Compute and store various metrics for the current network.
    final_results[name]['mean_segment_length'] = skeleton_paths['distance'].mean()
    final_results[name]['segment_length_list'] = skeleton_paths['distance'].values
    final_results[name]['mean_radius_length'] = skeleton_paths['radius_1'].mean()
    final_results[name]['radius_list'] = skeleton_paths['radius_1'].values
    final_results[name]['mean_segment_over_bounding_box_average_length'] = final_results[name]['mean_segment_length'] / bounding_box_volume**(1/3)

    # Determine cube side and starting proportion based on the network name.
    cube_side = shortest_side * 0.32
    starting_proportion = 5 if name != 'fruit_fly_3' else 4
    cube_side_list = [(starting_proportion / 100) * cube_side]

    # Determine the resolution for segments based on their size.
    segment_resolution = min(skeleton_paths['distance'].mean(), skeleton_paths['radius_1'].mean(), skeleton_paths['radius_2'].mean()) / 2
    final_results[name]['cube_side'] = cube_side_list[0]
    final_results[name]['segment_resolution_distance'] = segment_resolution

    print('Skeleton multiplication...')
    multiplied_skeleton_part_merged = []
    # Simplified skeleton multiplication process.
    for row in skeleton_paths[['pt_id_1', 'x_1', 'y_1', 'z_1', 'radius_1', 'pt_id_2', 'x_2', 'y_2', 'z_2', 'radius_2']].values:
        new_row = [list(row[1:5]) + list(row[6:10])]
        multiplied_skeleton_part_merged += new_row

    multiplied_skeleton = pd.DataFrame(multiplied_skeleton_part_merged, columns=['x_1', 'y_1', 'z_1', 'radius_1', 'x_2', 'y_2', 'z_2', 'radius_2'])
    skeleton_paths = multiplied_skeleton.copy()
    print('Skeleton multiplication complete.')

    # Compute the point cloud for skeleton paths.
    point_cloud_list = []
    for i, row in enumerate(skeleton_paths[['x_1', 'y_1', 'z_1', 'radius_1', 'x_2', 'y_2', 'z_2', 'radius_2']].values):
        print(f'**** Network: {name} | Computing fractal dimension {i} of {len(skeleton_paths)}')
        point_cloud_list += create_thick_line_point_list(row[:3], row[4:7], segment_resolution, row[3], row[7])
    point_cloud_skeleton_paths = pd.DataFrame(point_cloud_list, columns=['x', 'y', 'z'])
    
    # Load the paths data for linking paths and skeleton points.
    link_paths = pd.read_csv(path + name + '.paths.csv', index_col=[0])
    # Filter columns to include only the start and end points and path identifiers.
    link_paths = link_paths[['x_1', 'y_1', 'z_1', 'x_2', 'y_2', 'z_2', 'path_id']].copy()

    # Initialize lists to store computed densities.
    all_densities = []
    all_link_densities = []
    
    # Iterate over each cube side value to calculate densities within cubes of that edge length.
    for cube_side in cube_side_list:
        print('**** Network:', name, 'Cube side', cube_side, 'out of', cube_side_list)

        # Generate ranges for the 3D grid based on bounding box limits and cube side length.
        x_axis = np.arange(bounds[0], bounds[1] + cube_side, cube_side)
        y_axis = np.arange(bounds[2], bounds[3] + cube_side, cube_side)
        z_axis = np.arange(bounds[4], bounds[5] + cube_side, cube_side)

        # Initialize counters for filled and empty boxes in the grid.
        filled_boxes_count = 0
        empty_boxes_count = 0

        # Loop through each cube in the 3D grid.
        for i in range(len(x_axis) - 1):
            for j in range(len(y_axis) - 1):
                for k in range(len(z_axis) - 1):
                    # Identify points within the current cube.
                    points_inside = point_cloud_skeleton_paths[
                        (point_cloud_skeleton_paths['x'] >= x_axis[i]) &
                        (point_cloud_skeleton_paths['x'] < x_axis[i + 1]) &
                        (point_cloud_skeleton_paths['y'] >= y_axis[j]) &
                        (point_cloud_skeleton_paths['y'] < y_axis[j + 1]) &
                        (point_cloud_skeleton_paths['z'] >= z_axis[k]) &
                        (point_cloud_skeleton_paths['z'] < z_axis[k + 1])
                    ]

                    # Update counts based on whether the current cube contains points.
                    if len(points_inside) == 0:
                        empty_boxes_count += 1
                    else:
                        filled_boxes_count += 1
                        # Only compute densities for cubes with enough points to form a convex hull.
                        if len(points_inside) >= 4:
                            # Compute link density based on unique path IDs within the cube.
                            links_inside = link_paths[
                                (link_paths['x_1'] >= x_axis[i]) & (link_paths['x_1'] < x_axis[i + 1]) &
                                (link_paths['y_1'] >= y_axis[j]) & (link_paths['y_1'] < y_axis[j + 1]) &
                                (link_paths['z_1'] >= z_axis[k]) & (link_paths['z_2'] < z_axis[k + 1])
                            ]
                            link_density = 0 if len(links_inside) == 0 else len(set(links_inside['path_id'].values))
                            all_link_densities.append(link_density)

                            # Compute physical density using the volume of the convex hull divided by the cube volume.
                            convex_hull = ConvexHull(points_inside[['x', 'y', 'z']])
                            local_density = convex_hull.volume / cube_side**3
                            all_densities.append(local_density)
                            if local_density > 1:
                                print("Density larger than 1")

    # Store computed density values in the results dictionary.
    final_results[name]['density_list'] = all_densities
    final_results[name]['link_density_list'] = all_link_densities

    # Save the density results to a pickle file.
    with open(name + "_density_results.pkl", "wb") as h:
        pickle.dump(final_results[name], h)