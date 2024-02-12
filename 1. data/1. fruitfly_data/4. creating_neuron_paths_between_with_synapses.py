import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from neuprint import fetch_synapses, NeuronCriteria as NC, SynapseCriteria as SC
from neuprint import Client
from neuprint import fetch_adjacencies, NeuronCriteria as NC
from neuprint import fetch_neurons
from neuprint import merge_neuron_properties
import math
import json
from neuprint import fetch_synapse_connections
from neuprint import Client
from urllib.error import URLError, HTTPError
from urllib.request import urlopen
import timeit
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
import pickle

import numpy as np
import pandas as pd
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import pickle
import timeit
import itertools
from neuprint import fetch_roi_hierarchy


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from neuprint import fetch_synapses, NeuronCriteria as NC, SynapseCriteria as SC
from neuprint import Client
from neuprint import fetch_adjacencies, NeuronCriteria as NC
from neuprint import fetch_neurons
from neuprint import merge_neuron_properties
import math
import json
from neuprint import fetch_synapse_connections
from neuprint import Client
from urllib.error import URLError, HTTPError
from urllib.request import urlopen
import timeit
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
import pickle

import numpy as np
import pandas as pd
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import pickle
import timeit
import itertools
from neuprint import fetch_roi_hierarchy
import pandas as pd
from scipy.spatial import distance
import numpy as np
import networkx as nx

def add_bodyid(row,max_points,column):
    if len(row[column]) <= len(str(max_points)):
        return str(row['bodyId']) + '_' + str(row[column]) 
    else:
        return str(row[column])
        
def add_bodyid_bridges(row,max_points,column):
    if column =='pt_id1':
        if len(row[column]) <= len(str(max_points)):
            return str(row['bodyId_pre']) + '_' + str(row[column]) 
        else:
            return str(row[column])
    else:
        if len(row[column]) <= len(str(max_points)):
            return str(row['bodyId_post']) + '_' + str(row[column]) 
        else:
            return str(row[column])

        



TOKEN = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6Imx1a2FibGFnb2pldmljMTk5NUBnbWFpbC5jb20iLCJsZXZlbCI6Im5vYXV0aCIsImltYWdlLXVybCI6Imh0dHBzOi8vbGgzLmdvb2dsZXVzZXJjb250ZW50LmNvbS9hLS9BT2gxNEdqdDZpdFFGR2xTSWZUTElNTjRmcEt1QzZ3QmE2Rlp0WU1XYmpKV1ZBPXM5Ni1jP3N6PTUwP3N6PTUwIiwiZXhwIjoxODA5MTE1OTEyfQ.0h6CJp8xfQEpkW8a2_gqJUBrEA5GyBiZkNvDjRpoXoY" # <--- Paste your token here
           # (or define NEUPRINT_APPLICATION CREDENTIALS in your environment)
c = Client('neuprint.janelia.org', 'hemibrain:v1.2.1', TOKEN)
index_list = [1,2,3,4]
dataset_lengths = {}
neuron_df_dict = {}
    
for index_name in index_list:
    print('Dataset:',index_name)
    path_source_1 = '2. healed_neuron_skeletons/'
    path_source_2 = '3. synapse_connectome_information/'
    path_save = '4. points_segments_synapse_connected/'
    skeleton_points = pd.read_csv(path_source_1 + "fruit_fly_" + str(index_name) +'.points.csv')
    skeleton_segments = pd.read_csv(path_source_1 + "fruit_fly_" + str(index_name) +'.segments.csv')
    connectome_df =  pd.read_csv(path_source_2 + "fruit_fly_" + str(index_name) +'.connectome.csv',index_col=[0])
    #soma_locations_df = pd.read_csv(path_source_2 + "fruit_fly_" + str(index_name) +'_soma_locations.csv',index_col=[0])
    synapses_df = pd.read_csv(path_source_2 + "fruit_fly_" + str(index_name) +'_synapse_info.csv',index_col=[0])

    connectome_df = connectome_df[['bodyId_pre','bodyId_post']].drop_duplicates().copy()

    print("Adding synapse info")
    synapse_confidence = 0.9
    new_synapse_points_row_list = []
    new_synapse_segments_row_list = []
    new_synapse_bridge_row_list = []
    missing_synapse_pairs_list = []
    multiple_synapse_pairs_list = []
        
    max_points =   1
    count_missing = 0
    count_present = 0
    new_seg_id = 1
    distance_list = []
    for i,pair in enumerate(connectome_df.values):
        print(i,len(connectome_df.values))
        bodyId_pre = pair[0]
        bodyId_post = pair[1]
        row_df = synapses_df[(synapses_df['bodyId_pre'] == bodyId_pre) & (synapses_df['bodyId_post'] == bodyId_post) ].drop_duplicates()
        for row_index,row in row_df.iterrows():
                count_present += 1
                single_points_pre_df = skeleton_points[skeleton_points['bodyId'] == bodyId_pre].copy()
                single_points_post_df = skeleton_points[skeleton_points['bodyId'] == bodyId_post].copy()
                single_points_pre_coordinates_df = single_points_pre_df[['x','y','z']].values
                single_points_post_coordinates_df = single_points_post_df[['x','y','z']].values
                radius_pre = single_points_pre_df['radius'].min()
                radius_post = single_points_post_df['radius'].min()
                new_point_pre =  max_points  + 1
                max_points =  max_points  + 1
                new_point_post = max_points + 1
                max_points  =  max_points + 1
                point_pre_coordinates = list(row[['x_pre','y_pre','z_pre']].values)
                point_post_coordinates = list(row[['x_post','y_post','z_post']].values)
                
                #Adding synapse points
                new_row_pre = pd.DataFrame({'pt_id':[new_point_pre],'x':[point_pre_coordinates[0]],'y':[point_pre_coordinates[1]],'z':[point_pre_coordinates[2]],'radius':[radius_pre],'diameter':[radius_pre*2],'bodyId':[bodyId_pre]})
                new_row_post = pd.DataFrame({'pt_id':[new_point_post],'x':[point_post_coordinates[0]],'y':[point_post_coordinates[1]],'z':[point_post_coordinates[2]],'radius':[radius_post],'diameter':[radius_post*2],'bodyId':[bodyId_post]})
                new_synapse_points_row_list.append( new_row_pre)
                new_synapse_points_row_list.append( new_row_post)
                distance_list.append(distance.euclidean(point_pre_coordinates ,point_post_coordinates))
                #Adding pre synapse segment
                dist_matrix = distance.cdist([point_pre_coordinates],single_points_pre_coordinates_df,'euclidean')
                distance_list.append(dist_matrix .min())
                closest_point_index = np.array(np.where(dist_matrix  == dist_matrix .min())).flatten()[1]
                row_pre = single_points_pre_df.iloc[closest_point_index]
                closest_point_pt_id  = row_pre['pt_id']
                new_seg_id_pre = new_seg_id
                new_seg_id += 1
                new_segment_pre = pd.DataFrame({'seg_id':[new_seg_id_pre],'pt_id1':[new_point_pre],'pt_id2':[ closest_point_pt_id ],'bodyId':[bodyId_pre]})
                new_synapse_segments_row = new_segment_pre.copy()
                
                
                #Adding post synapse segment
                dist_matrix = distance.cdist([point_post_coordinates],single_points_post_coordinates_df,'euclidean')
                distance_list.append(dist_matrix .min())
                closest_point_index = np.array(np.where(dist_matrix  == dist_matrix .min())).flatten()[1]
                row_post = single_points_post_df.iloc[closest_point_index]
                closest_point_pt_id  = row_post['pt_id']
                new_seg_id_post =  new_seg_id
                new_seg_id += 1
                new_segment_post = pd.DataFrame({'seg_id':[new_seg_id_post],'pt_id1':[new_point_post],'pt_id2':[ closest_point_pt_id ],'bodyId':[bodyId_post]})
                new_synapse_segments_row = pd.concat([new_synapse_segments_row,new_segment_post ])
                new_synapse_segments_row_list.append(new_synapse_segments_row)
                
                #Creating bridge segment
                new_segment_bridge = pd.DataFrame({'seg_id':[new_seg_id],'pt_id1':[new_point_pre],'pt_id2':[new_point_post],'bodyId_pre':[bodyId_pre],'bodyId_post':[bodyId_post]})
                new_seg_id += 1
                new_synapse_bridge_row = new_segment_bridge.copy()
                new_synapse_bridge_row_list.append(new_synapse_bridge_row)
                

    points_to_add = pd.concat(new_synapse_points_row_list)
    segments_to_add = pd.concat(new_synapse_segments_row_list)
    bridges_to_add = pd.concat(new_synapse_bridge_row_list)


    print("Adding new bodyids")
    save_points = skeleton_points.copy()
    save_segments = skeleton_segments.copy()

    #skeleton_segments['pt_id2'] = skeleton_segments['pt_id2'].copy().astype('int64')
    skeleton_segments['pt_id1_bodyid'] = skeleton_segments['pt_id1'].copy()
    skeleton_segments['pt_id2_bodyid'] =  skeleton_segments['pt_id2'].copy()
    skeleton_segments['bodyId'] = skeleton_segments['bodyId'].copy().astype('int64')
    skeleton_segments['bodyId_pre'] = skeleton_segments['bodyId'].copy().astype('int64')
    skeleton_segments['bodyId_post'] = skeleton_segments['bodyId'].copy().astype('int64')
    skeleton_points['pt_id_bodyid'] = skeleton_points['pt_id'].copy() 
    
    #Points are added in reverse
    points_to_add['pt_id'] = points_to_add['pt_id'].copy().astype('int64')
    points_to_add['bodyId'] = points_to_add['bodyId'].copy().astype('int64')
    points_to_add['pt_id_bodyid'] = points_to_add['bodyId'].astype(str).copy()  + '_' +   points_to_add['pt_id'].astype(str).copy()
    skeleton_points = skeleton_points.append(points_to_add)
    
    segments_to_add['pt_id1'] = segments_to_add['pt_id1'].copy().astype(str)
    segments_to_add['pt_id2'] = segments_to_add['pt_id2'].copy().astype(str)
    segments_to_add['bodyId'] = segments_to_add['bodyId'].copy().astype('int64')
    segments_to_add['bodyId_pre'] = segments_to_add['bodyId'].copy().astype('int64')
    segments_to_add['bodyId_post'] = segments_to_add['bodyId'].copy().astype('int64')
    segments_to_add['pt_id2_bodyid'] = segments_to_add.apply(lambda x:add_bodyid(x,max_points,'pt_id2'),axis=1)
    segments_to_add['pt_id1_bodyid'] = segments_to_add.apply(lambda x:add_bodyid(x,max_points,'pt_id1'),axis=1)
    skeleton_segments = skeleton_segments.append(segments_to_add)
    
    bridges_to_add['pt_id1'] = bridges_to_add['pt_id1'].copy().astype(str)
    bridges_to_add['pt_id2'] =bridges_to_add['pt_id2'].copy().astype(str)
    bridges_to_add['bodyId_pre'] =bridges_to_add['bodyId_pre'].copy().astype('int64')
    bridges_to_add['bodyId_post'] =bridges_to_add['bodyId_post'].copy().astype('int64')
    bridges_to_add['pt_id1_bodyid'] = bridges_to_add.apply(lambda x:add_bodyid_bridges(x,max_points,'pt_id1'),axis=1)
    bridges_to_add['pt_id2_bodyid'] =bridges_to_add.apply(lambda x:add_bodyid_bridges(x,max_points,'pt_id2'),axis=1)
    skeleton_segments = skeleton_segments.append(bridges_to_add)
    
    skeleton_points_df = skeleton_points[['pt_id_bodyid','x','y','z','diameter','bodyId']].copy()
    skeleton_points_df.rename(columns={'pt_id_bodyid':'pt_id','bodyId':'neuron_id'},inplace=True)
    skeleton_points_dict = {}
    for i,row in enumerate(skeleton_points_df.values):
        skeleton_points_dict[row[0]] = i
    skeleton_points_df['pt_id'] = skeleton_points_df['pt_id'].apply(lambda x:skeleton_points_dict[x]).copy()
    skeleton_points_df.to_csv(path_save + 'fruit_fly_'+str(index_name) +'.points.csv')
    
    skeleton_segments['seg_id'] = np.arange(0,len(skeleton_segments))
    skeleton_segments_df = skeleton_segments[['seg_id','pt_id1_bodyid','pt_id2_bodyid','bodyId_pre','bodyId_post']].copy()
    skeleton_segments_df.rename(columns={'pt_id1_bodyid':'pt_id1','pt_id2_bodyid':'pt_id2','bodyId_pre':'neuron_id_pre','bodyId_post':'neuron_id_post'},inplace=True)
    skeleton_segments_df.set_index('seg_id', inplace = True)
    skeleton_segments_df['pt_id1'] = skeleton_segments_df['pt_id1'].apply(lambda x:skeleton_points_dict[x]).copy()
    skeleton_segments_df['pt_id2'] = skeleton_segments_df['pt_id2'].apply(lambda x:skeleton_points_dict[x]).copy()
    skeleton_segments_df.drop_duplicates(inplace=True)
    skeleton_segments_df.to_csv(path_save + 'fruit_fly_'+str(index_name) +'.segments.csv')
    print('Maximum distance of new segments',np.max(distance_list))
    print('Maximum degree aded to a single point',segments_to_add['pt_id2'].value_counts()[0])