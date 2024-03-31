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

from physical_properties_functions import *


# Loop through each network name and perform computations
name_list = ['human_neuron','rat_neuron',
           'monkey_neuron','zebrafish_neuron', 'vascular_2','vascular_3','vascular_1','mitochondrial','anthill','root_1','root_2','fruit_fly_2','fruit_fly_3','fruit_fly_1','fruit_fly_4']
final_results = {}

for name in name_list:
    print('**** Network:', name)
    print('Starting computation')
    final_results[name] = {}
    path = '../1. data/3. final_data/'
    
    # Load and preprocess skeleton paths data
    skeleton_paths = pd.read_csv(path + name + '.paths.csv', index_col=[0])
    skeleton_paths.drop_duplicates(inplace=True)

    # Compute distances and volumes for each path
    calc_distances_and_volumes(skeleton_paths)

    # Store some calculated metrics
    final_results[name]['median_segment_length'] = skeleton_paths['distance'].median()
    final_results[name]['median_segment_radius'] = np.median((skeleton_paths['radius_1'].values + skeleton_paths['radius_2'].values)/2)
    final_results[name]['radius_list'] = (skeleton_paths['radius_1'].values + skeleton_paths['radius_2'].values)/2
    final_results[name]['segment_length_list'] = skeleton_paths['distance'].values
    
    # Compute complementary straightness and update results
    final_results[name].update(compute_complementary_straightness(skeleton_paths))

    # Save the updated results to a pickle file
    with open(f"0. link_properties/{name}_basic_results.pkl", "wb") as h:
        pickle.dump(final_results[name], h)