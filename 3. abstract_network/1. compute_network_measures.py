import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
from networkx.algorithms.connectivity import local_edge_connectivity
import networkx as nx
import seaborn as sns    

# Define a list of network names for analysis.
name_list = ['human_neuron', 'rat_neuron', 'monkey_neuron', 'zebrafish_neuron', 'vascular_2', 'vascular_3', 'vascular_1', 'mitochondrial', 'anthill', 'root_1', 'root_2', 'fruit_fly_2', 'fruit_fly_3', 'fruit_fly_1', 'fruit_fly_4']

# Loop through each network name to perform analyses.
for name in name_list:
    print('******', name)
    network_measures_dict = {}  # Initialize a dictionary to store network measures.
    path_source = '../1. data/3. final_data/'  # Define the base path for data files.
    
    # Load paths data from CSV file for the current network.
    link_paths = pd.read_csv(path_source + name + '.paths.csv', index_col=[0])
    
    # Extract unique path, source, and target information.
    path_bodyid_list = link_paths[['path_id', 'source', 'target']].drop_duplicates().values.tolist()
    
    G = nx.MultiGraph()  # Initialize a MultiGraph object to represent the network.
    path_bodyid_dict = {}  # Initialize a dictionary to map path IDs to source-target pairs.
    bodyid_path_dict = {}  # Initialize a dictionary to map source-target pairs to path IDs.

    # Populate dictionaries with path information.
    for path_bodyid in path_bodyid_list:
        path_bodyid_dict[path_bodyid[0]] = (path_bodyid[1], path_bodyid[2])
        
        # Update mappings in both directions for source-target pairs and path IDs.
        if (path_bodyid[1], path_bodyid[2]) not in bodyid_path_dict:
            bodyid_path_dict[(path_bodyid[1], path_bodyid[2])] = [path_bodyid[0]]
            bodyid_path_dict[(path_bodyid[2], path_bodyid[1])] = [path_bodyid[0]]
        else:
            bodyid_path_dict[(path_bodyid[1], path_bodyid[2])].append(path_bodyid[0])
            bodyid_path_dict[(path_bodyid[2], path_bodyid[1])].append(path_bodyid[0])
    
    # Add edges to the graph based on the source-target pairs.
    for path_id, bodyid_edge in path_bodyid_dict.items():
        G.add_edge(bodyid_edge[0], bodyid_edge[1])

    print('Calculating degrees')
    degrees = dict(nx.degree(G))  # Calculate the degree of each node.

    # Initialize a dictionary to store the sum of degrees for each path.
    link_degree = {}
    for path_id in path_bodyid_dict.keys():
        node_pair = path_bodyid_dict[path_id]
        link_degree[path_id] = degrees[node_pair[0]] + degrees[node_pair[1]]  # Sum the degrees of source and target nodes for each path.

    print('Calculating network diameter')
    network_diameter = nx.diameter(G)  # Calculate the diameter of the network.

    print('Edge betweenness centrality')
    betw_dict = dict(nx.edge_betweenness_centrality(G))  # Calculate edge betweenness centrality for all edges.

    # Initialize a dictionary to map betweenness centrality values to path IDs.
    betw_path_dict = {}
    for path_id in path_bodyid_dict.keys():
        node_pair = path_bodyid_dict[path_id]
        # Handle the case where the node pair might be reversed in the centrality dictionary.
        if node_pair in betw_dict:
            betw_path_dict[path_id] = betw_dict[node_pair]
        else:
            betw_path_dict[path_id] = betw_dict[(node_pair[1], node_pair[0])]

    # Store computed measures in the dictionary.
    network_measures_dict['network_diameter'] = network_diameter
    network_measures_dict['link_degree_dict'] = link_degree
    network_measures_dict['betw_dict'] = betw_path_dict

    # Save the network measures dictionary to a pickle file for later use.
    with open("1. network_measures_results/" + name + "_network_measures_dict.pkl", "wb") as h:
        pickle.dump(network_measures_dict, h)