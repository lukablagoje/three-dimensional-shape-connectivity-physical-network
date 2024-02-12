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
import networkx as nx
from scipy.spatial import distance

TOKEN = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6Imx1a2FibGFnb2pldmljMTk5NUBnbWFpbC5jb20iLCJsZXZlbCI6Im5vYXV0aCIsImltYWdlLXVybCI6Imh0dHBzOi8vbGgzLmdvb2dsZXVzZXJjb250ZW50LmNvbS9hLS9BT2gxNEdqdDZpdFFGR2xTSWZUTElNTjRmcEt1QzZ3QmE2Rlp0WU1XYmpKV1ZBPXM5Ni1jP3N6PTUwP3N6PTUwIiwiZXhwIjoxODA5MTE1OTEyfQ.0h6CJp8xfQEpkW8a2_gqJUBrEA5GyBiZkNvDjRpoXoY" # <--- Paste your token here
           # (or define NEUPRINT_APPLICATION CREDENTIALS in your environment)

c = Client('neuprint.janelia.org', 'hemibrain:v1.2.1', TOKEN)
#skeletons_dict = {1:'POC_96', 4: 'GF(R)_10',3:'mALT(L)_20',2:'GC_28'}
skeletons_dict = {4:'GF(R)_10'}
path_source = '1. neuron_regions_skeleton/'
path_save = '3. synapse_connectome_information/'
all_body_dict = {}
for i, skeleton_region in skeletons_dict.items():
    print('Looking for synapses')
    name = 'fruit_fly_'+str(i)
                   
    region_df = pd.read_csv(path_source + skeleton_region + '.csv')
    bodyid_list= list(set(region_df['bodyId'].values))
    print('Looking for synapses')
    neuron_df, connectome_df = fetch_adjacencies(bodyid_list,bodyid_list)
    connectome_df.to_csv(path_save + name +'.connectome.csv')
    neuron_pairs = connectome_df[['bodyId_pre','bodyId_post']].values
    synapse_df_list = []
    count = 0
    count_true = 0
    for pair in neuron_pairs:
        print(count_true,len(neuron_pairs))
        try:
            eb_conns = fetch_synapse_connections( NC(bodyId=[pair[0]]),NC(bodyId=[pair[1]]), None)
            synapse_df_list.append(eb_conns)
        except:
            try:
                eb_conns = fetch_synapse_connections( NC(bodyId=[pair[1]]),NC(bodyId=[pair[0]]), None)
                synapse_df_list.append(eb_conns)
            except:
                print("No synapse pair found")
        count +=1
        count_true += 1
        if count % 100: 
            print("Completed",count,len(neuron_pairs))
            synapse_df = pd.concat(synapse_df_list)
            synapse_df.to_csv(path_save + name + '_synapse_info.csv')
            count = 0
    synapse_df.to_csv(path_save + name + '_synapse_info.csv')