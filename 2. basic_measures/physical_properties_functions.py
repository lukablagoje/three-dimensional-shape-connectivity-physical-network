import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ast
import math
from scipy.spatial import distance as dst
from scipy.spatial import ConvexHull
from scipy.stats import linregress


def create_line_point_radii_list(start_point, end_point, radius_start, radius_end, density, no_radius):
    """
    Create a list of points along a line with varying radii.

    Parameters:
    :param start_point: The starting point of the line (a list or array of x, y, z coordinates).
    :param end_point: The ending point of the line (a list or array of x, y, z coordinates).
    :param radius_start: The radius at the starting point of the line.
    :param radius_end: The radius at the ending point of the line.
    :param density: The number of points to generate along the line.
    :param no_radius: A boolean indicating whether to include radius in the output (True to exclude radius).

    Returns:
    :return: A list of points along the line. Each point is a list of x, y, z coordinates, 
             optionally including the radius if no_radius is False.
    """

    # Generate a list of parameters for interpolating along the line
    parameter_list = np.linspace(0, 1, density)
    line_point_list = []

    # Iterating over each parameter to generate points along the line
    for t in parameter_list:
        # Interpolate the x, y, and z coordinates
        x = (end_point[0] - start_point[0]) * t + start_point[0]
        y = (end_point[1] - start_point[1]) * t + start_point[1]
        z = (end_point[2] - start_point[2]) * t + start_point[2]
        
        # Interpolate the radius
        radius = (radius_end - radius_start) * t + radius_start

        # Add the point to the list, with or without radius based on no_radius flag
        if no_radius:
            line_point_list.append([x, y, z])
        else:
            line_point_list.append([x, y, z, radius])

    return line_point_list


def create_inner_circles(center_point, radius, distance_between_circles, x, y):
    """
    Create inner circles for larger radii.
    
    :param center_point: Center point of the circle.
    :param radius: Radius of the circle.
    :param distance_between_circles: Distance between circles.
    :param x: Orthogonal unit vector to the line.
    :param y: Another orthogonal unit vector to the line.
    :return: List of points forming inner circles.
    """
    inner_points = []

    # Check if the radius is large enough to require inner circles
    if 2 * radius > distance_between_circles:
        # Calculate the number of inner circles based on radius and distance between circles
        number_of_inner_circles = max(math.ceil(2 * radius / distance_between_circles), 3)
        inner_radius_list = [radius * parameter for parameter in np.linspace(0, 1, number_of_inner_circles)[1:-1]]

        # Generate points for each inner circle
        for inner_radius in inner_radius_list:
            n_vertices = max(math.ceil(2 * inner_radius * np.pi / distance_between_circles), 3)
            for j in range(n_vertices):
                # Calculate the coordinates of each point on the inner circle
                angle = (j / n_vertices) * 2 * np.pi
                circle_x = center_point[0] + inner_radius * np.cos(angle) * x[0] + inner_radius * np.sin(angle) * y[0]
                circle_y = center_point[1] + inner_radius * np.cos(angle) * x[1] + inner_radius * np.sin(angle) * y[1]
                circle_z = center_point[2] + inner_radius * np.cos(angle) * x[2] + inner_radius * np.sin(angle) * y[2]
                inner_points.append([circle_x, circle_y, circle_z])

    return inner_points

def handle_additional_points(start_point, orientation_point, radius_1, radius_2, distance_between_circles, x, y):
    """
    Handle additional points for larger radii and create inner circles if necessary.
    
    :param start_point: Starting point of the line.
    :param orientation_point: Orientation point of the line.
    :param radius_1: Radius at the start point.
    :param radius_2: Radius at the orientation point.
    :param distance_between_circles: Distance between circles along the line.
    :param x: Orthogonal unit vector to the line.
    :param y: Another orthogonal unit vector to the line.
    :return: List of additional points.
    """
    additional_points = []

    # Check if radii are larger than the distance between circles.
    if radius_1 > distance_between_circles:
        additional_points.append(start_point)
    if radius_2 > distance_between_circles:
        additional_points.append(orientation_point)

    # Create inner circles for larger radii at the start and orientation points.
    additional_points += create_inner_circles(start_point, radius_1, distance_between_circles, x, y)
    additional_points += create_inner_circles(orientation_point, radius_2, distance_between_circles, x, y)

    return additional_points


def create_thick_line_point_list(start_point, orientation_point, distance_between_circles, radius_1, radius_2):
    """
    Generates a list of points representing a thick line between two points in 3D space, with variable radii.

    This function creates a "thick" line (cylinder or tapered cylinder) by generating points in circles perpendicular to the line defined by start_point and orientation_point. The circles' radii change linearly from radius_1 to radius_2 along the length of the line. This method is useful for visualizing tubes or pipes in 3D space with varying thickness.

    Parameters:
    - start_point (array-like): The starting point of the line segment in 3D space (x, y, z).
    - orientation_point (array-like): The ending point of the line segment in 3D space (x, y, z), which defines the direction of the line.
    - distance_between_circles (float): The distance between each generated circle along the line. Smaller values result in a smoother appearance.
    - radius_1 (float): The radius of the circle at the starting point. This can be used to create tapered effects.
    - radius_2 (float): The radius of the circle at the orientation point, allowing for tapered effects.
    
    Returns:
    - list: A list of 3D points (x, y, z) that, when connected, form a thick line or tube between the start and orientation points.
    
    The function calculates the direction vector between the start and orientation points, normalizes it, and then generates two orthogonal vectors to form the basis for the circular sections. Points on these circles are calculated using trigonometric functions and added to the output list. Additional handling is performed to ensure that the tube has a consistent appearance, even when the radii are significantly different.
    """
    # Calculate the Euclidean distance between the start and orientation points.
    distance = np.linalg.norm(np.array(orientation_point) - np.array(start_point))
    # Determine the number of points to create along the line, ensuring at least 2 points.
    number_of_points = max(int(distance / distance_between_circles), 2)

    # Create a line point list using radii information.
    no_radius =True
    line_point_list = create_line_point_radii_list(start_point, orientation_point, radius_1, radius_2,number_of_points,no_radius)
    thick_points_list = []

    # Calculate unit vector along the line.
    k = np.array(orientation_point) - np.array(start_point)
    k /= np.linalg.norm(k)

    # Generate two orthogonal unit vectors to 'k'.
    x = np.random.randn(3)
    x -= x.dot(k) * k
    x /= np.linalg.norm(x)
    y = np.cross(k, x)

    # Iterate over each point in the line and create circles perpendicular to the line.
    for i, center_point in enumerate(line_point_list):
        parameter = i / number_of_points
        radius = (radius_2 - radius_1) * parameter + radius_1
        n_vertices = max(math.ceil(2 * radius * np.pi / distance_between_circles), 3)
        for j in range(n_vertices):
            angle = (j / n_vertices) * 2 * np.pi
            circle_x = center_point[0] + radius * np.cos(angle) * x[0] + radius * np.sin(angle) * y[0]
            circle_y = center_point[1] + radius * np.cos(angle) * x[1] + radius * np.sin(angle) * y[1]
            circle_z = center_point[2] + radius * np.cos(angle) * x[2] + radius * np.sin(angle) * y[2]
            v = [circle_x, circle_y, circle_z]
            thick_points_list.append(v)

    # Handle additional points for larger radii and create inner circles if needed.
    thick_points_list += handle_additional_points(start_point, orientation_point, radius_1, radius_2, distance_between_circles, x, y)

    return thick_points_list + line_point_list
    

def bounding_box_from_merged_skeleton(merged_skeleton_dataset):
    """
    Calculate the bounding box of a merged skeleton dataset.

    The bounding box is defined by the minimum and maximum coordinates 
    in each of the x, y, and z dimensions.

    :param merged_skeleton_dataset: A dataset containing skeleton data with x, y, z coordinates.
    :return: A list containing the minimum and maximum coordinates in the format [x_min, x_max, y_min, y_max, z_min, z_max].
    """

    # Extract x, y, and z coordinates from the dataset
    x_coords = merged_skeleton_dataset[['x_1', 'x_2']].values
    y_coords = merged_skeleton_dataset[['y_1', 'y_2']].values
    z_coords = merged_skeleton_dataset[['z_1', 'z_2']].values

    # Calculate the minimum and maximum values for each dimension
    x_min, x_max = np.min(x_coords), np.max(x_coords)
    y_min, y_max = np.min(y_coords), np.max(y_coords)
    z_min, z_max = np.min(z_coords), np.max(z_coords)

    return [x_min, x_max, y_min, y_max, z_min, z_max]

def perform_skeleton_multiplication(skeleton_paths, cube_side):
    """
    Perform skeleton multiplication to generate an expanded set of skeleton parts.

    Parameters:
    :param skeleton_paths: DataFrame containing skeleton path data.
    :param cube_side: The length of the side of the cube used for skeleton multiplication.

    Returns:
    :return: A list containing the expanded set of skeleton parts.
    """
    multiplied_skeleton_part_merged = []

    for row in skeleton_paths[['pt_id_1', 'x_1', 'y_1', 'z_1', 'radius_1', 'pt_id_2', 'x_2', 'y_2', 'z_2', 'radius_2']].values:
        start_point = row[1:4]
        radius_start = row[4]
        end_point = row[6:9]
        radius_end = row[9]
        point_to_point_dist = dst.euclidean(start_point, end_point)
        ratio = (point_to_point_dist + radius_start + radius_end) / cube_side

        if ratio > 0.2:
            multiplication_factor = 7  # Previously was 2 + 5
            new_points = create_line_point_radii_list(start_point, end_point, radius_start, radius_end, multiplication_factor, False)
            for i in range(multiplication_factor - 1):
                multiplied_skeleton_part_merged += [list(new_points[i]) + list(new_points[i + 1])]
        else:
            new_row = [list(start_point) + [radius_start] + list(end_point) + [radius_end]]
            multiplied_skeleton_part_merged += new_row

    return multiplied_skeleton_part_merged

def generate_point_cloud_list(skeleton_paths, segment_resolution):
    """
    Generate a point cloud list from the given skeleton paths.

    Parameters:
    :param skeleton_paths: DataFrame containing the multiplied skeleton path data.
    :param segment_resolution: The resolution of the segments in the skeleton paths.

    Returns:
    :return: A list representing the point cloud generated from the skeleton paths.
    """
    point_cloud_list = []

    for row in skeleton_paths[['x_1', 'y_1', 'z_1', 'radius_1', 'x_2', 'y_2', 'z_2', 'radius_2']].values:
        start_point = row[0:3]
        radius_start = row[3]
        end_point = row[4:7]
        radius_end = row[7]
        point_cloud_list += create_thick_line_point_list(start_point, end_point, segment_resolution, radius_start, radius_end)

    return point_cloud_list

def process_network_data(name, path):
    """
    Process the network data for a given name and path.

    :param name: Name of the network.
    :param path: Path to the network data.
    :return: Dictionary containing various calculated metrics for the network.
    """
    print('**** Network:', name)
    print('Starting computation')
    network_results = {}

    # Load and preprocess skeleton paths data
    skeleton_paths = pd.read_csv(f'{path}{name}.paths.csv', index_col=0)
    skeleton_paths.drop_duplicates(inplace=True)

    # Calculate bounding box and related metrics
    bounds = bounding_box_from_merged_skeleton(skeleton_paths)
    bounding_box_volume = np.prod(np.diff(bounds.reshape(3, 2), axis=1).flatten())
    longest_side = np.max(np.diff(bounds.reshape(3, 2), axis=1))
    shortest_side = np.min(np.diff(bounds.reshape(3, 2), axis=1))

    # Compute distances and volumes for each path
    skeleton_paths['distance'] = np.linalg.norm(skeleton_paths[['x_2', 'y_2', 'z_2']].values - skeleton_paths[['x_1', 'y_1', 'z_1']].values, axis=1)
    skeleton_paths['volume'] = (1/3) * np.pi * np.prod(skeleton_paths[['radius_1', 'radius_2', 'radius_1']].values, axis=1) * skeleton_paths['distance']

    # Store some calculated metrics
    network_results['mean_segment_length'] = skeleton_paths['distance'].mean()
    network_results['segment_length_list'] = skeleton_paths['distance'].values
    network_results['mean_radius_length'] = skeleton_paths['radius_1'].mean()
    network_results['radius_list'] =(skeleton_paths['radius_1'].values + skeleton_paths['radius_2'].values)/2
    network_results['mean_segment_over_bounding_box_average_length'] = network_results['mean_segment_length'] / (bounding_box_volume ** (1/3))

    # Define cube side and segment resolution
    cube_side = shortest_side * 0.32
    cube_side_list = [proportion / 100 * cube_side for proportion in np.arange(15, 105, 5)]
    segment_resolution = min(skeleton_paths['distance'].mean(), skeleton_paths['radius_1'].mean(), skeleton_paths['radius_2'].mean()) / 2
    network_results['cube_side'] = cube_side
    network_results['cube_side_list'] = cube_side_list
    network_results['segment_resolution_distance'] = segment_resolution

def calc_distances_and_volumes(skeleton_paths):
    """
    Calculate distances and volumes for each path segment in the skeleton paths.
    :param skeleton_paths: DataFrame containing skeleton path data.
    """
    skeleton_paths['distance'] = np.sqrt(np.sum(np.square(skeleton_paths[['x_2', 'y_2', 'z_2']].values - skeleton_paths[['x_1', 'y_1', 'z_1']].values), axis=1))
    skeleton_paths['volume'] = (1/3) * np.pi * (skeleton_paths['radius_1']**2 + skeleton_paths['radius_2']**2 + skeleton_paths['radius_1']*skeleton_paths['radius_2'])*skeleton_paths['distance']

def compute_complementary_straightness(skeleton_paths):
    """
    Compute the complementary straightness for each path in the skeleton paths.
    :param skeleton_paths: DataFrame containing skeleton path data.
    :return: Dictionary containing calculated straightness and related metrics.
    """
    results = {'c_straightness': {}, 'link_volume': {}, 'link_segments': {}, 'link_path_length': {},'link_aspect_ratio':{},'link_straightness':{}}
    
    for path_id in skeleton_paths['path_id'].unique():
        path_segments = skeleton_paths[skeleton_paths['path_id'] == path_id]
        total_path_length = path_segments['distance'].sum()
        total_volume = path_segments['volume'].sum()
        mean_link_radius = np.mean([path_segments.iloc[0][['radius_1']].values] + list(path_segments[['radius_2']].values))
        
        start_loc = path_segments.iloc[0][['x_1', 'y_1', 'z_1']].values
        end_loc = path_segments.iloc[-1][['x_2', 'y_2', 'z_2']].values

        euclidean_distance = np.linalg.norm(end_loc - start_loc)
        straightness = euclidean_distance / total_path_length if total_path_length > 0 else 0
        results['link_aspect_ratio'][path_id] = mean_link_radius/total_path_length
        results['link_straightness'][path_id] = straightness
        results['c_straightness'][path_id] = 1 - straightness
        results['link_volume'][path_id] = total_volume
        results['link_segments'][path_id] = len(path_segments)
        results['link_path_length'][path_id] = total_path_length
        
    total_volume_of_network = sum(results['link_volume'].values())
    results['link_volume_normed'] = {path_id: vol / total_volume_of_network for path_id, vol in results['link_volume'].items()}

    return results