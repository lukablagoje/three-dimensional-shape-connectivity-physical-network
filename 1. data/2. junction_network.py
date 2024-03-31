import pandas as pd
import networkx as nx

def label_paths(segments):
    """
    Labels each segment of a skeleton with a unique path identifier based on connectivity.
    
    Args:
    segments (np.array): An array of segment connections, where each row contains the IDs
                         of the connected points (pt_id1 and pt_id2).
    
    Returns:
    tuple: Contains three elements:
        - A dictionary mapping each edge to a path identifier.
        - A dictionary mapping each path identifier to its start and end node.
        - A dictionary mapping each path identifier to a list of its edges.
    """
    G = nx.Graph()
    G.add_edges_from(segments)  # Create a graph from the segments

    # Check if the graph is connected
    if not nx.is_connected(G):
        print("Graph not connected")
        return 0

    # Identify junction nodes: nodes with a degree not equal to 2
    junctions = [node for node in G.nodes() if G.degree(node) != 2]

    path_id = 0
    path_labels = {}
    path_endpoints = {}  # Dictionary to store start and end nodes of each path
    path_edgelists = {}
    visited_edges = set()  # Set to keep track of visited edges

    for junction in junctions:
        # Explore each neighboring node of the junction
        for neighbor in G.neighbors(junction):
            edge = tuple(sorted((junction, neighbor)))
            if edge not in visited_edges:
                # Start a new path from the junction
                start_node = junction
                current_path = [edge]
                visited_edges.add(edge)
                next_node = neighbor

                # Extend the path through nodes of degree 2
                while G.degree(next_node) == 2:
                    next_neighbors = list(G.neighbors(next_node))
                    for next_neighbor in next_neighbors:
                        next_edge = tuple(sorted((next_node, next_neighbor)))
                        if next_edge not in visited_edges:
                            current_path.append(next_edge)
                            visited_edges.add(next_edge)
                            next_node = next_neighbor
                            break

                end_node = next_node  # End of the path
                # Store the edges and endpoints of the path
                path_edgelists[path_id] = current_path
                path_endpoints[path_id] = (start_node, end_node)
                # Label each edge with the current path ID
                for edge in current_path:
                    path_labels[edge] = path_id
                path_id += 1

    return path_labels, path_endpoints, path_edgelists

def sort_paths(path_edgelists, path_endpoints):
    """
    Sorts the edges within each path to ensure they are in order from start to end node.
    
    Args:
    path_edgelists (dict): A dictionary where each key is a path ID and each value is a list of edges in that path.
    path_endpoints (dict): A dictionary mapping each path ID to its start and end node.
    
    Returns:
    dict: A dictionary where each key is a path ID and each value is a list of sorted edges from start to end node.
    """
    sorted_paths = {}

    for path_id, edges in path_edgelists.items():
        start_node, end_node = path_endpoints[path_id]
        G = nx.Graph()
        G.add_edges_from(edges)
        # Compute the shortest path from start to end node within the path's graph
        path = nx.shortest_path(G, source=start_node, target=end_node)
        # Convert the node path back into a list of edges
        path_edges = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
        sorted_paths[path_id] = path_edges

    return sorted_paths

name_list = ['human_neuron','rat_neuron',
            'monkey_neuron','zebrafish_neuron',
    'vascular_2','vascular_3','vascular_1','mitochondrial','root_1','root_2','anthill']
# Loop through each structure in the provided list
for name in name_list:
    print('*******', name)
    
    # Define source and destination paths for reading and saving files
    path_source = '3. healed_skeletons/'
    path_save = '4. intersection_connectomes/'
    
    # Load the segments data for the current structure
    segments_df = pd.read_csv(path_source + name + '.segments.csv', index_col=0).astype(int)
    
    # Label paths in the segment data, which includes identifying path start and end points, and edge lists for each path
    path_labels, path_endpoints, path_edgelists = label_paths(segments_df[['pt_id1', 'pt_id2']].astype(int).values)
    
    # Sort the paths to ensure they are ordered from start node to end node
    sorted_path_edgelists = sort_paths(path_edgelists, path_endpoints)
    
    # Initialize a dictionary to hold the combined segment and point data
    df_dict = {'pt_id_1': [], 'pt_id_2': [], 'path_id': [], 'source': [], 'target': []}
    
    # Loop through each path to populate the dictionary with segment and path endpoint data
    for path_id in sorted_path_edgelists:
        for edge in sorted_path_edgelists[path_id]:
            if edge[0] != edge[1]:  # Ensure the edge is not a loop back to the same point
                df_dict['pt_id_1'].append(edge[0])
                df_dict['pt_id_2'].append(edge[1])
                df_dict['path_id'].append(path_id)
                df_dict['source'].append(path_endpoints[path_id][0])  # Add the source node of the path
                df_dict['target'].append(path_endpoints[path_id][1])  # Add the target node of the path
    
    # Convert the combined segment and path data into a DataFrame
    df_segment_labels = pd.DataFrame(df_dict)

    # Load the points data for the current structure
    points_df = pd.read_csv(path_source + name + '.points.csv')
    points_df['radius'] = points_df['diameter'] / 2  # Convert diameter to radius
    points_df = points_df[['pt_id', 'x', 'y', 'z', 'radius']]  # Select relevant columns

    # Reinitialize the dictionary to combine segment, point, and path information
    df_dict = {'pt_id_1': [], 'x_1': [], 'y_1': [], 'z_1': [], 'radius_1': [], 'pt_id_2': [], 'x_2': [], 'y_2': [], 'z_2': [], 'radius_2': [], 'path_id': [], 'source': [], 'target': []}
    
    # Loop through each row in the segments DataFrame to populate the dictionary with corresponding point data
    for i, row in df_segment_labels.iterrows():
        if i % 1000 == 0:
            print('Combining segment, point and path information.', i, len(df_segment_labels))
        
        # Append data for the first point of the segment
        df_dict['pt_id_1'].append(row['pt_id_1'])
        point_source_df = points_df[points_df['pt_id'] == row['pt_id_1']]
        df_dict['x_1'].append(point_source_df['x'].values[0])
        df_dict['y_1'].append(point_source_df['y'].values[0])
        df_dict['z_1'].append(point_source_df['z'].values[0])
        df_dict['radius_1'].append(point_source_df['radius'].values[0])

        # Append data for the second point of the segment
        df_dict['pt_id_2'].append(row['pt_id_2'])
        point_target_df = points_df[points_df['pt_id'] == row['pt_id_2']]
        df_dict['x_2'].append(point_target_df['x'].values[0])
        df_dict['y_2'].append(point_target_df['y'].values[0])
        df_dict['z_2'].append(point_target_df['z'].values[0])
        df_dict['radius_2'].append(point_target_df['radius'].values[0])

        # Append path identifier and source-target information
        df_dict['path_id'].append(row['path_id'])
        df_dict['source'].append(row['source'])
        df_dict['target'].append(row['target'])

    # Save the combined data to a CSV file for further use
    pd.DataFrame(df_dict).to_csv(path_save + name + '.paths.csv')