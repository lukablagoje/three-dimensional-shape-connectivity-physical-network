import numpy as np
import pandas as pd
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import pickle
import timeit
import itertools
from neuprint import fetch_roi_hierarchy
import networkx as nx
from scipy.spatial import distance

skeletons_dict = {1:'POC_96', 4: 'GF(R)_10',3:'mALT(L)_20',2:'GC_28'}
path_source = '1. neuron_regions_skeleton/'
path_save = ''#'2. healed_neuron_skeletons/'
for df_index,skeleton_region in skeletons_dict.items():
    #skeleton_region = 'ME(R)_3721'
    df = pd.read_csv(path_source + skeleton_region + '.csv',index_col=[0])
    name = skeleton_region
    #df = pd.read_csv(path_source+ name_original,delimiter=' ')
    segments = df[['link','rowId','bodyId']].copy()
    segments = segments.rename(columns = {'link':'pt_id1', 'rowId':'pt_id2'})
    segments = segments[(segments['pt_id1'] != -1) & (segments['pt_id2'] != -1)]
    segments.reset_index(inplace=True)
    segments = segments.rename(columns = {'level_0':'seg_id'})
    segments.drop(['index'],inplace=True,axis=1)
    #segments.set_index(['seg_id'],inplace=True)
    #segments.to_csv(path + name + '.segments.csv')

    points = df[['rowId','x','y','z','radius','bodyId']].copy()
    points = points.rename(columns = {'rowId':'pt_id'})
    points['diameter'] = points['radius'] * 2
    points.reset_index()
    #points.to_csv(path + name+'.points.csv')

    bodyid_components_dict = {}
    healed_list_segments = []
    healed_list_points = []
    for i,bodyid in enumerate(list(set(segments['bodyId'].values))):
        print('Neuron number',i,'out of',len(list(set(segments['bodyId'].values))))
        segments_df = segments[segments['bodyId'] == bodyid].copy()
        points_df = points[points['bodyId'] == bodyid].copy()
        disconnected_skeleton = True
        iteration_step = 1
        while disconnected_skeleton:
            print('Iteration step',iteration_step)
            #Create a graph from segments and check if it's fully connected
            G = nx.Graph()
            G.add_edges_from(segments_df[['pt_id1','pt_id2']].values)
            components = list(nx.connected_components(G))
            print('The number of disconnected components is',len(components))
            if len(components) == 1:
                disconnected_skeleton = False
                break
            #Creating a dictionary for points and their coordinates and radius
            node_coordinate_dict = {}
            node_coordinate_radius_dict = {}
            component_dict = {}
            component_coordinate_dict = {}
            for i,component in enumerate(list(nx.connected_components(G))):
                component_dict[i] = component
                for node in component:
                    coordinates = list(points_df[points_df['pt_id'] == node][['x','y','z']].values)
                    node_coordinate_dict[node] = list(coordinates[0])
                    node_coordinate_radius_dict[node] = list(points_df[points_df['pt_id'] == node][['x','y','z','radius']].values[0])
                component_coordinate_dict[i] = [node_coordinate_dict[node] for node in component]

            #Finding distance between components
            component_closest_dict = {}
            for component_1 in component_dict.keys():
                for component_2 in component_dict.keys():
                    #print(component_1,component_2)
                    if component_1 != component_2 and not((component_1,component_2) in component_closest_dict.keys()) and not((component_2,component_1) in component_closest_dict.keys()):
                        component_closest_dict[(component_1,component_2)]= np.amin(distance.cdist(component_coordinate_dict[component_1],  component_coordinate_dict[component_2], 'euclidean'))

            #Finding closests components
            target_components = {}
            for component_1 in component_dict.keys():
                minimum_value = np.max(list(component_closest_dict.values())) +1
                for component_2 in component_dict.keys():
                    if (component_1,component_2) in list(component_closest_dict.keys()) or (component_2,component_1) in list(component_closest_dict.keys()):
                        if (component_1,component_2) in list(component_closest_dict.keys()):
                            current_value = component_closest_dict[(component_1,component_2)] 
                        if (component_2,component_1) in list(component_closest_dict.keys()):
                            current_value = component_closest_dict[(component_2,component_1)] 
                        if current_value < minimum_value:
                            minimum_value = current_value
                            current_component = component_2
                if not(current_component in target_components.keys()):
                    target_components[(component_1,current_component)] = minimum_value

            #Finding minimum distance nodes to connect
            nodes_to_be_healed = []
            for pair in target_components.keys():
                #print(pair,len(list(target_components.keys())))
                node_list_1 = component_dict[pair[0]]
                node_list_2 = component_dict[pair[1]]
                minimum_distance = np.max(list(component_closest_dict.values())) +1
                for node_1 in node_list_1:
                    node_1_coordinates = node_coordinate_dict[node_1] 
                    for node_2 in node_list_2:
                        node_2_coordinates = node_coordinate_dict[node_2] 
                        node_distance = distance.euclidean(node_1_coordinates,node_2_coordinates)
                        if node_distance < minimum_distance:
                            minimum_distance = node_distance
                            minimum_distance_nodes = (node_1,node_2)
                nodes_to_be_healed.append(minimum_distance_nodes)

            # Creatings segment between targeted nodes
            for node_pair in nodes_to_be_healed:
                node_1 = node_pair[0]
                node_2 = node_pair[1]
                row = [node_1] + list(node_coordinate_radius_dict[node_1]) + [node_2] + list(node_coordinate_radius_dict[node_2])
                row = pd.DataFrame({'pt_id1':[node_1],'pt_id2':[node_2]},index=[np.max(segments_df.index) +1])
                segments_df = pd.concat([segments_df, row], axis=0)
        print('Healing complete')
        segments_df['bodyId'] = bodyid
        segments_df['seg_id'] = np.arange(0,len(segments_df))
        segments_df.set_index('seg_id', inplace = True)
        segments_df.drop_duplicates(inplace=True)
        segments_df['pt_id1'] = segments_df['pt_id1'].astype(int).astype(str) + '_' + segments_df['bodyId'].astype(str)
        segments_df['pt_id2'] = segments_df['pt_id2'].astype(int).astype(str) + '_' + segments_df['bodyId'].astype(str)
        points_df['pt_id'] = points_df['pt_id'].astype(int).astype(str) + '_' + points_df['bodyId'].astype(str)
        points_df.set_index('pt_id',inplace=True)
        healed_list_segments.append(segments_df)
        healed_list_points.append(points_df)
    full_skeleton_segments = pd.concat(healed_list_segments)
    full_skeleton_segments.reset_index(inplace=True)
    full_skeleton_segments.to_csv(path_save +'fruit_fly_' + str(df_index)+ '.segments.csv')
    full_skeleton_points = pd.concat(healed_list_points)
    full_skeleton_points.to_csv(path_save +'fruit_fly_' + str(df_index)+ '.points.csv')