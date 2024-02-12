import pandas as pd
import networkx as nx
def label_paths(segments):
    G = nx.Graph()
    G.add_edges_from(segments)
    # Check if the graph is connected
    if not  nx.is_connected(G):
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
        for neighbor in G.neighbors(junction):
            edge = tuple(sorted((junction, neighbor)))
            if edge not in visited_edges:
                start_node = junction  # Starting node of the path
                current_path = [edge]
                visited_edges.add(edge)
                next_node = neighbor

                while G.degree(next_node) == 2:
                    next_neighbors = list(G.neighbors(next_node))
                    for next_neighbor in next_neighbors:
                        next_edge = tuple(sorted((next_node, next_neighbor)))
                        if next_edge not in visited_edges:
                            current_path.append(next_edge)
                            visited_edges.add(next_edge)
                            next_node = next_neighbor
                            break

                end_node = next_node  # Ending node of the path
                path_edgelists[path_id] = []
                for edge in current_path:
                    path_labels[edge] = path_id
                    path_edgelists[path_id].append(edge)
                path_endpoints[path_id] = (start_node, end_node)  # Store the start and end nodes
                path_id += 1

    return path_labels, path_endpoints,path_edgelists


def sort_paths(path_edgelists, path_endpoints):
    sorted_paths = {}

    for path_id, edges in path_edgelists.items():
        if path_id % 1000 == 0: 
            print('Sorting path',path_id,len(path_edgelists))
        start_node, end_node = path_endpoints[path_id]
        G = nx.Graph()
        # Initialize with the first edge that contains the start node
        G.add_edges_from(path_edgelists[path_id])
        path = nx.shortest_path(G, source=start_node, target=end_node)
        path_edges = [(path[i], path[i + 1]) for i in range(len(path) - 1)]

        sorted_paths[path_id] = path_edges 
    return sorted_paths


name_list = ['human_neuron','rat_neuron',
            'monkey_neuron','zebrafish_neuron',
    'vascular_2','vascular_3','vascular_1','mitochondrial','root_1','root_2','anthill']
for name in name_list:
    print('*******',name)
    path_source = '3. healed_skeletons/'
    path_save = '4. intersection_connectomes/'
    segments_df = pd.read_csv(path_source + name+'.segments.csv',index_col=[0]).astype(int)
    segments_df
    path_labels, path_endpoints, path_edgelists = label_paths(segments_df[['pt_id1','pt_id2']].astype(int).values)
    sorted_path_edgelists = sort_paths(path_edgelists, path_endpoints)
    df_dict = {'pt_id_1': [], 'pt_id_2': [], 'path_id': [], 'source': [], 'target': []}
    for path_id in sorted_path_edgelists:
        for edge in sorted_path_edgelists[path_id]:
            # Remove self-loops
            if edge[0] != edge[1]:
                df_dict['pt_id_1'].append(edge[0])
                df_dict['pt_id_2'].append(edge[1])
                df_dict['path_id'].append(path_id)
                df_dict['source'].append( path_endpoints[path_id][0])
                df_dict['target'].append( path_endpoints[path_id][1])
    df_segment_labels = pd.DataFrame(df_dict)

    points_df = pd.read_csv(path_source + name + '.points.csv')
    points_df['radius'] = points_df['diameter']/2
    points_df = points_df[['pt_id','x','y','z','radius']]

    df_dict = {'pt_id_1': [], 'x_1':[],'y_1':[],'z_1':[],'radius_1':[], 'pt_id_2': [], 'x_2':[],'y_2':[],'z_2':[],'radius_2':[],'path_id': [], 'source': [], 'target': []}
    for i, row in df_segment_labels.iterrows():
        if i%1000 == 0:
            print('Combining segment, point and path information.',i,len(df_segment_labels))
        df_dict['pt_id_1'].append(row['pt_id_1'])
        point_source_df = points_df[points_df['pt_id'] == row['pt_id_1']]
        df_dict['x_1'].append(point_source_df['x'].values[0])
        df_dict['y_1'].append(point_source_df['y'].values[0])
        df_dict['z_1'].append(point_source_df['z'].values[0])
        df_dict['radius_1'].append(point_source_df['radius'].values[0])

        df_dict['pt_id_2'].append(row['pt_id_2'])
        point_target_df = points_df[points_df['pt_id'] == row['pt_id_2']]
        df_dict['x_2'].append(point_target_df['x'].values[0])
        df_dict['y_2'].append(point_target_df['y'].values[0])
        df_dict['z_2'].append(point_target_df['z'].values[0])
        df_dict['radius_2'].append(point_target_df['radius'].values[0])

        df_dict['path_id'].append(row['path_id'])
        df_dict['source'].append(row['source'])
        df_dict['target'].append(row['target'])

    pd.DataFrame(df_dict).to_csv(path_save + name + '.paths.csv')