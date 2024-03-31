import numpy as np
import pandas as pd
import pickle

def e_distance(start_point, end_point):
    """
    Calculates the Euclidean distance between two points in 3D space.
    
    Args:
    start_point (list): The starting point coordinates [x, y, z].
    end_point (list): The ending point coordinates [x, y, z].
    
    Returns:
    float: The Euclidean distance between the start and end points.
    """
    distance = np.sqrt((end_point[0] - start_point[0])**2 + (end_point[1] - start_point[1])**2 + (end_point[2] - start_point[2])**2)
    return distance

def rw_length(rw):
    """
    Computes the total length of a random walk (rw) by summing the Euclidean distances between consecutive points.
    
    Args:
    rw (list): A list of points representing the random walk, where each point is [x, y, z].
    
    Returns:
    float: The total length of the random walk, rounded to 4 decimal places.
    """
    total_distance = 0
    for i in range(len(rw)-1):
        start = rw[i][0:3]
        end = rw[i+1][0:3]
        total_distance += np.linalg.norm(np.array(end) - np.array(start))
        
    return np.round(total_distance, 4)


# List of biological structure names to process
name_list = ['human_neuron', 'rat_neuron', 'monkey_neuron', 'zebrafish_neuron', 'vascular_2', 'vascular_3', 'vascular_1', 'mitochondrial', 'root_1', 'root_2', 'anthill']

# Directories for source data and where to save the refined paths
path_source = '4. intersection_connectomes/'
path_save = '5. intersection_connectomes_skeleton_resolution/'

# Dictionary to store clustering steps information for each structure
clustering_steps_dict = {}
for name in name_list:
    print('*****', name)
    # Load the paths CSV as a DataFrame
    df = pd.read_csv(path_source + name + '.paths.csv', index_col=[0])
    
    # Initialize dictionaries to store new path information and clustering stopping steps
    new_path_dict = {}
    stopping_step_dict = {}

    # Extract unique path IDs from the DataFrame
    path_id_list = list(set(df['path_id'].values))
    
    # Iterate over each path ID to refine path segment connections
    for count, path_id in enumerate(path_id_list):
        print(count, len(path_id_list))

        # Extract path data for the current path ID
        df_one_path = df[df['path_id'] == path_id]
        bodyid_pre = df_one_path['source'].values[0]  # Starting node ID of the path
        bodyid_post = df_one_path['target'].values[0]  # Ending node ID of the path

        rw = []  # Initialize a list to store path data including points and radii

        # Populate rw with point data from the DataFrame
        for point in df_one_path[['x_1', 'y_1', 'z_1', 'pt_id_1', 'radius_1']].values:
            rw.append(np.array(point))
        rw.append(df_one_path[['x_2', 'y_2', 'z_2', 'pt_id_2', 'radius_2']].values[-1])  # Add the last point separately

        # Calculate the total length of the path
        total_path_length = rw_length(rw)

        # Initialize lists to store metrics used for evaluating the clustering process
        mean_cos_angle_list = []
        std_cos_angle_list = []
        removed_cos_angle_list = []
        segment_list = []  # Initialize a list to store segments and their directional vectors

        # Iterate through each segment in the path to calculate its direction vector
        for index, segment in enumerate(rw):
            if index < len(rw) - 1:
                start = rw[index][:3]
                end = rw[index + 1][:3]
                diff_vector = np.array(end) - np.array(start)  # Calculate the difference vector

                if np.linalg.norm(diff_vector) != 0:  # Exclude zero-length segments
                    diff_vector_hat = diff_vector / np.linalg.norm(diff_vector)  # Normalize the difference vector
                    segment_list.append(list(rw[index]) + list(rw[index + 1]) + list(diff_vector_hat))

        # Calculate the number of clustering steps as one less than the number of segments
        clustering_steps = len(segment_list) - 2

        # Start the clustering process to refine path segment connections
        stopping_criteria_search = True
        # Iterate through each clustering step to refine path connections
        for step in range(0, clustering_steps + 2):
            print('Clustering step', step, 'out of', clustering_steps + 1)

            # Reinitialize the segment list for the current clustering step
            segment_list = []

            # Iterate through each segment in the path
            for index, segment in enumerate(rw):
                # Ensure we're not at the last segment to calculate directional vectors
                if index < len(rw) - 1:
                    # Extract start and end points of the current segment
                    start = rw[index][:3]
                    end = rw[index + 1][:3]
                    # Calculate the difference vector and normalize it to get the direction vector
                    diff_vector = np.array(end) - np.array(start)
                    diff_vector_hat = diff_vector / np.linalg.norm(diff_vector)
                    # Append the segment information along with its direction vector to the list
                    segment_list.append(rw[index].tolist() + rw[index + 1].tolist() + diff_vector_hat.tolist())

            # Convert the segment list into a DataFrame for easier manipulation
            df_segments = pd.DataFrame(segment_list, columns=['x_1', 'y_1', 'z_1', 'pt_id_1', 'radius_1', 'x_2', 'y_2', 'z_2', 'pt_id_2', 'radius_2', 'x_norm', 'y_norm', 'z_norm'])

            # Initialize a dictionary to hold the cosine angle between consecutive direction vectors
            cos_angle_dict = {}

            # Calculate cosine angles for direction vectors
            for i, norm_diff_vec in enumerate(df_segments[['x_norm', 'y_norm', 'z_norm']].values):
                if i < len(df_segments) - 1:
                    # Cosine angle calculation between two consecutive direction vectors
                    cos_angle = np.dot(norm_diff_vec, df_segments.iloc[i + 1][['x_norm', 'y_norm', 'z_norm']].values)
                    # Apply special treatment for the first and last segments to avoid altering their positions
                    if i == 0 or i == len(df_segments) - 2:
                        cos_angle_dict[(i, i + 1)] = 0  # Assign zero to avoid clustering these segments
                    else:
                        cos_angle_dict[(i, i + 1)] = cos_angle  # Assign the calculated cosine angle

            # Sort segments based on their cosine angle in descending order
            sorted_segment_indices = sorted(cos_angle_dict, key=cos_angle_dict.get, reverse=True)

            # Begin segment clustering based on the sorted cosine angles
            if len(cos_angle_dict) > 0:
                # Record the mean and standard deviation of cosine angles before removing the segment with the largest angle
                mean_cos_angle_list.append(np.mean(list(cos_angle_dict.values())))
                std_cos_angle_list.append(np.std(list(cos_angle_dict.values())))
                removed_cos_angle_list.append(cos_angle_dict[sorted_segment_indices[0]])

                # Identify and merge the segments with the maximum cosine angle
                index_1, index_2 = sorted_segment_indices[0]
                first_row = df_segments.iloc[index_1]
                second_row = df_segments.iloc[index_2]
                new_start = first_row[['x_1', 'y_1', 'z_1', 'pt_id_1', 'radius_1']].values
                new_end = second_row[['x_2', 'y_2', 'z_2', 'pt_id_2', 'radius_2']].values
                new_diff_vec = new_end[:3] - new_start[:3]
                new_diff_vec_hat = new_diff_vec / np.linalg.norm(new_diff_vec)
                # Drop the original segments and insert the merged segment
                df_dropped = df_segments.drop(index=[index_1, index_2])
                df_dropped.loc[(index_1 + index_2) / 2] = list(new_start) + list(new_end) + list(new_diff_vec_hat)
            else:
                df_dropped = df_segments.copy()

            # Reset the DataFrame index after dropping segments
            df_reset = df_dropped.sort_index().reset_index(drop=True)

            # Prepare the path data for the next iteration based on updated segments
            rw_new = [np.array(point) for point in df_reset[['x_1', 'y_1', 'z_1', 'pt_id_1', 'radius_1']].values]
            rw_new.append(df_reset[['x_2', 'y_2', 'z_2', 'pt_id_2', 'radius_2']].values[-1])

            # Calculate the length ratio of the updated path to the original path
            length_ratio = rw_length(rw_new) / total_path_length
            length_ratio_list.append(length_ratio)

            # Check if the clustering process significantly altered the path length
            if length_ratio < 0.9999999 and stopping_criteria_search:
                # Record the step at which significant alteration occurred and stop the search
                stopping_step = step
                stopping_criteria_search = False
                stopping_step_dict[path_id] = step
                break

            # Update the path data for the next clustering step
            rw = rw_new

        # If the search wasn't stopped explicitly, record the final step
        if stopping_criteria_search:
            stopping_step_dict[path_id] = step

        # Save the refined path for later use
        new_path_dict[(path_id, bodyid_pre, bodyid_post)] = rw
    # Store the stopping step information for the current structure
    clustering_steps_dict[name] = stopping_step_dict
    
    # Compile refined path information into a new DataFrame
    row_list = []
    for path_key, points in new_path_dict.items():
        for i in range(len(points)-1):
            # Extract segment and point data for each connection within the path
            row = list(points[i][:3]) + [str(int(points[i][3]))] + [points[i][4]] + list(points[i+1][:3]) + [str(int(points[i+1][3]))] + [points[i+1][4]] + list(path_key)
            row_list.append(row)
    
    new_df = pd.DataFrame(row_list, columns=['x_1', 'y_1', 'z_1', 'pt_id_1', 'radius_1', 'x_2', 'y_2', 'z_2', 'pt_id_2', 'radius_2', 'path_id', 'source', 'target'])
    
    # Save the refined paths to a new CSV file
    new_df.to_csv(path_save + name + '.paths.csv')
