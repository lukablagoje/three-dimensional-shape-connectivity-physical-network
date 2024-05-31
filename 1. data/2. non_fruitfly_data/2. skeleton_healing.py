import pandas as pd
import numpy as np
import networkx as nx
from scipy.spatial.distance import cdist
from itertools import combinations


def heal_skeletons(segments, points):
    # Convert points data to a dictionary of node attributes
    node_attributes = points.set_index('pt_id')[['x', 'y', 'z', 'radius']].to_dict('index')

    # Initialize the graph
    G = nx.Graph()
    G.add_edges_from(segments[['pt_id1', 'pt_id2']].values)

    # Ensure seg_id is an integer and find the maximum seg_id for initializing new ids
    segments['seg_id'] = segments['seg_id'].astype(int)
    max_seg_id = segments['seg_id'].max()

    # Initialize the iteration count
    iteration_number = 0

    # Continue healing until all components are connected
    while True:
        print('Iteration number', iteration_number + 1)
        components = list(nx.connected_components(G))
        if len(components) == 1:
            print("Skeleton is fully connected.")
            break

        print("Number of components", len(components))
        healed = False  # Flag to check if any healing was done in the loop

        # Calculate distances for the first pair of components that can be connected
        for component1, component2 in combinations(components, 2):
            component1_coords = np.array([[node_attributes[node]['x'], node_attributes[node]['y'], node_attributes[node]['z']] for node in component1])
            component2_coords = np.array([[node_attributes[node]['x'], node_attributes[node]['y'], node_attributes[node]['z']] for node in component2])
            distances = cdist(component1_coords, component2_coords)
            min_dist_idx = np.unravel_index(np.argmin(distances), distances.shape)
            closest_nodes = (list(component1)[min_dist_idx[0]], list(component2)[min_dist_idx[1]])

            # Create and add the new healing segment
            new_row = pd.DataFrame({'seg_id': [max_seg_id + 1],
                                    'pt_id1': [closest_nodes[0]],
                                    'pt_id2': [closest_nodes[1]],
                                    'segment_healing_done': ['healed']})
            segments = pd.concat([segments, new_row], ignore_index=True)
            G.add_edge(closest_nodes[0], closest_nodes[1])
            max_seg_id += 1  # Increment the maximum segment ID
            healed = True

            # Break after adding one edge
            break

        iteration_number += 1
        if not healed:
            print("No more healing possible.")
            break

    print("Healing complete. Number of components:", len(list(nx.connected_components(G))))
    return segments, points

name_list = ['human_neuron', 'rat_neuron', 'monkey_neuron', 'zebrafish_neuron', 'vascular_2', 'vascular_3', 'vascular_1', 'mitochondrial', 'root_1', 'root_2', 'anthill']
path_source = '2. points_segments/'
path_save = '3. healed_skeletons/'
for name in name_list:
    print('*****', name)
    points = pd.read_csv(path_source + name +'.points.csv')
    segments = pd.read_csv(path_source + name+'.segments.csv')
    segments['segment_healing_done'] = 'original'
    points['radius'] = points['diameter']/2
    segments_healed, points_healed = heal_skeletons(segments, points)
    segments_healed.to_csv(path_save + name + '.segments.csv')
    points_healed.to_csv(path_save + name + '.points.csv')