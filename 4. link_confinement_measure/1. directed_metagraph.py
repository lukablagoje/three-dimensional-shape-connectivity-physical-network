import pandas as pd
#from vedo import Line, show,merge,Points
import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
import math
import os
import timeit
import random
def create_line_point_list(start_point, end_point, density):

    parameter_list = np.linspace(0,1,density)
    line_point_list = []

    for t in parameter_list:
            x = (end_point[0] - start_point[0]) * t + start_point[0]
            y = (end_point[1] - start_point[1]) * t + start_point[1]
            z = (end_point[2] - start_point[2]) * t + start_point[2]
            line_point_list.append([x, y, z])

    return line_point_list


def create_thick_line_point_list(start_point,orientation_point,distance_between_circles,radius_1,radius_2):
    distance = np.sqrt((orientation_point[0]-start_point[0])**2 +(orientation_point[1]-start_point[1])**2 +(orientation_point[2]-start_point[2])**2)
    number_of_points =  np.max([int(distance/distance_between_circles),2])
    #number_of_points = density
    if radius_1 == 0 and radius_2 == 0:
        line_point_list = create_line_point_list(start_point,orientation_point,number_of_points)
        return line_point_list 
    else:
        line_point_list = create_line_point_list(start_point,orientation_point,number_of_points)
        thick_points_list = []
        k = np.array([orientation_point[0]-start_point[0], orientation_point[1]-start_point[1] ,orientation_point[2]-start_point[2] ],dtype = float )
        k /= np.linalg.norm(k)
        x = np.random.randn(3)  # take a random vector
        x -= x.dot(k) * k       # make it orthogonal to k
        x /= np.linalg.norm(x)  # normalize it
        y = np.cross(k, x)
        for i, center_point in enumerate(line_point_list):
                parameter = i/len(line_point_list)
                radius = (radius_2 - radius_1) * parameter + radius_1
                n_vertices =  np.max([math.ceil(2 * radius * np.pi/ distance_between_circles),3])
                circle_point_list = []
                for j in range(0,n_vertices):
                    angle = (j / n_vertices) * 2 * np.pi

                    circle_x = center_point[0] + radius * np.cos(angle)*x[0] +  radius * np.sin(angle)*y[0]#* np.cos(incline_angle)
                    circle_y = center_point[1] + radius * np.cos(angle)*x[1] +  radius * np.sin(angle)*y[1]
                    circle_z = center_point[2] + radius * np.cos(angle)*x[2] +  radius * np.sin(angle)*y[2]

                    v = [circle_x,circle_y, circle_z]
                    #circle_point_list.append(v)
                    thick_points_list.append(v)

        if radius_1 > distance_between_circles:
            thick_points_list.append(start_point)
        if radius_2 > distance_between_circles:
            thick_points_list.append(orientation_point)

        if 2 * radius_1 > distance_between_circles:
            number_of_inner_circles = np.max([math.ceil( 2 * radius_1 / distance_between_circles),3])
            inner_radius_list = [radius_1 * parameter for parameter in np.linspace(0,1, number_of_inner_circles)[1:-1]]

            for inner_radius in inner_radius_list:
                center_point = start_point
                radius = inner_radius
                n_vertices =  np.max([math.ceil(2 * radius * np.pi/ segment_resolution),3])
                circle_point_list = []
                for j in range(0,n_vertices):
                    angle = (j / n_vertices) * 2 * np.pi

                    circle_x = center_point[0] + radius * np.cos(angle)*x[0] +  radius * np.sin(angle)*y[0]#* np.cos(incline_angle)
                    circle_y = center_point[1] + radius * np.cos(angle)*x[1] +  radius * np.sin(angle)*y[1]
                    circle_z = center_point[2] + radius * np.cos(angle)*x[2] +  radius * np.sin(angle)*y[2]

                    v = [circle_x,circle_y, circle_z]
                    #circle_point_list.append(v)
                    thick_points_list.append(v)
        if 2 * radius_2 > distance_between_circles:
            number_of_inner_circles = np.max([math.ceil( 2 * radius_2 / distance_between_circles),3])
            inner_radius_list = [radius_2 * parameter for parameter in np.linspace(0,1, number_of_inner_circles)[1:-1]]
            for inner_radius in inner_radius_list:
                center_point = orientation_point
                radius = inner_radius
                n_vertices =  np.max([math.ceil(2 * radius * np.pi/ segment_resolution),3])
                circle_point_list = []
                for j in range(0,n_vertices):
                    angle = (j / n_vertices) * 2 * np.pi

                    circle_x = center_point[0] + radius * np.cos(angle)*x[0] +  radius * np.sin(angle)*y[0]#* np.cos(incline_angle)
                    circle_y = center_point[1] + radius * np.cos(angle)*x[1] +  radius * np.sin(angle)*y[1]
                    circle_z = center_point[2] + radius * np.cos(angle)*x[2] +  radius * np.sin(angle)*y[2]

                    v = [circle_x,circle_y, circle_z]
                    #circle_point_list.append(v)
                    thick_points_list.append(v)

        return thick_points_list + line_point_list 




def randomize_segments_with_thickness(segment_list,radius_list,seed):
    diff_vector_list = []
    for segment in segment_list:
        diff_vector = np.array(segment[1]) - np.array(segment[0])
        diff_vector_list.append(diff_vector)
    starting_point = segment_list[0][0]
    np.random.seed(seed)   
    a = diff_vector_list
    b = radius_list
    c = list(zip(a, b))

    np.random.shuffle(c)
    
    diff_vector_list,radius_list = zip(*c)
    randomized_segment_list = []
    

    for diff_vector in diff_vector_list:
        new_second_point= starting_point + diff_vector
        randomized_segment_list.append([np.array(starting_point),np.array(new_second_point)])
        starting_point = new_second_point
    return randomized_segment_list, radius_list

def intersecting_bodyid_from_kd_tree(intersections,dict_of_all_points):
    point_indices = set( [num for sublist in intersections for num in sublist])
    all_points_keys = list( dict_of_all_points.keys())
    bodyid_intersected =  []
    for points_index in point_indices:
        for i in range(len(all_points_keys)):
            if all_points_keys[i] >= points_index:
                bodyid = all_points_keys[i]
                bodyid_intersected.append(dict_of_all_points[bodyid])
                break
    return bodyid_intersected
    
    
def list_dict_all_points(points_bodyid):

    list_of_all_points = []
    dict_of_all_points = {}
    for key in list(points_bodyid.keys()):

        list_of_all_points += points_bodyid[key]
        end = (len(list_of_all_points)) - 1
        dict_of_all_points[end] = key
    return list_of_all_points, dict_of_all_points
path_save = '3. directed_metagraph_results/'
#path_source = '../1. data/5. skeleton_resolution/'
path_source = '../1. data/3. final_data/'
#name_list = ['human_neuron','rat_neuron','monkey_neuron',
#'zebrafish_neuron','vascular_2','vascular_3','vascular_1','mitochondrial','anthill','root_1','root_2']
name_list = ['human_neuron','rat_neuron',
            'monkey_neuron','zebrafish_neuron', 'vascular_2','vascular_3','vascular_1','mitochondrial','root_1','root_2','anthill']
for name in name_list:
    print(name)
    skeleton_paths = pd.read_csv(path_source + name+'.paths.csv',index_col=[0])
    skeleton_paths['distance'] = np.sqrt((skeleton_paths['x_2']-  skeleton_paths['x_1'])**2 +  
                                               (skeleton_paths['y_2']-  skeleton_paths['y_1'])**2 +
                                               (skeleton_paths['z_2']-  skeleton_paths['z_1'])**2 )
    skeleton_paths.drop_duplicates(inplace=True)
    radius_mean_1 = skeleton_paths['radius_1'].mean()
    radius_mean_2 = skeleton_paths['radius_2'].mean()
    segment_distance_mean = skeleton_paths['distance'].mean()
    segment_resolution = np.min([segment_distance_mean,radius_mean_1,radius_mean_2 ])

    segments_path_dict = {}
    radius_path_dict = {}

    path_id_list = list(set(skeleton_paths['path_id'].values.tolist()))
    #path_id_list = path_id_list
    for i, path_id in enumerate(path_id_list):
        #print(i,len(path_id_list))
        segment_list = []
        for segment_joined in skeleton_paths[skeleton_paths['path_id'] == path_id][['x_1','y_1','z_1','x_2','y_2','z_2']].values:
             segment_list.append([ segment_joined[0:3],segment_joined[3:]  ]   )
        radius_list = []
        for radius_joined in skeleton_paths[skeleton_paths['path_id'] == path_id][['radius_1','radius_2']].values:
              radius_list.append([ radius_joined[0],radius_joined[1] ]   )
        segments_path_dict[path_id] =  segment_list
        radius_path_dict[path_id] = radius_list


    #path_id_list = list(set(skeleton_paths['path_id'].values.tolist()))
    path_bodyid_list = list(skeleton_paths[['path_id','bodyId_pre','bodyId_post']].values)


    kd_tree_original = {}
    for count,path_id in enumerate(path_id_list):
        print('Creating original KD-trees',count,len(path_id_list))
        segment_list = segments_path_dict[path_id]
        radius_list = radius_path_dict[path_id]
        thick_line_list = []
        for i in range(len(segment_list)):
             if not(np.array_equal(segment_list[i][0],segment_list[i][1])):
                 thick_line_list += create_thick_line_point_list(segment_list[i][0],segment_list[i][1],segment_resolution,radius_list[i][0],radius_list[i][1])
        kd_tree_original[path_id] = thick_line_list


        #print('Time for',count,'/',len(path_id_list) ,'is',endtime - startime)

    #with open('kd_trees_storage/'+"points_bodyid_distance_" + str(segment_resolution) + ".pkl", "wb") as h:
    #    pickle.dump(kd_tree_original, h)

    path_bodyid_list = skeleton_paths[['path_id','bodyId_pre','bodyId_post']].values.tolist()
    all_body_id_list = list(set(list(set(skeleton_paths['bodyId_pre'].values.tolist())) + list(set(skeleton_paths['bodyId_post'].values.tolist())) ))
    path_bodyid_dict = {}
    for bodyid in all_body_id_list:
        path_bodyid_dict[bodyid] = []
    for path_body in path_bodyid_list:
        path_bodyid_dict[path_body[1]].append(path_body[0])
        path_bodyid_dict[path_body[2]].append(path_body[0])
        
    for key in path_bodyid_dict.keys():
        path_bodyid_dict[key] = list(set(path_bodyid_dict[key]))


    adjacent_paths = {}
    for path_body in path_bodyid_list:
        path_id_1 = path_body[0]
        path_list_1 = path_bodyid_dict[path_body[1]] 
        path_list_2 = path_bodyid_dict[path_body[2]]
        path_id_list_adj = list(set(path_list_1 + path_list_2))
        path_id_list_adj.remove(path_id_1)
        adjacent_paths[path_id_1] = path_id_list_adj
        #paths_to_exclude_dict[path_body[0]].remove(path_body[0])
    #with open('kd_trees_storage/'+ "paths_to_exclude_dict.pkl", "wb") as h:
    #    pickle.dump(paths_to_exclude_dict, h)        

    print('T list start')

    t_list = [1]
    t_dict_results = {}
    N_trials = 20
    sample_size_paths = 500
    start_trial = 0
    len_list = len(path_id_list)
    for t in t_list:
        t_dict_results[t] = {}
        path_id_list_sample = path_id_list
        kd_tree_original = {}
        for count,path_id in enumerate(path_id_list):
            print(count,len(path_id_list))
            segment_list = segments_path_dict[path_id]
            radius_list = radius_path_dict[path_id]
            thick_line_list = []
            for i in range(len(segment_list)):
                 if not(np.array_equal(segment_list[i][0],segment_list[i][1])):
                     thick_line_list += create_thick_line_point_list(segment_list[i][0],segment_list[i][1],segment_resolution,radius_list[i][0]*t,radius_list[i][1]*t)
            kd_tree_original[path_id] = thick_line_list
        list_of_all_points, dict_of_all_points = list_dict_all_points(kd_tree_original)
        kd_tree_large = KDTree(list_of_all_points,balanced_tree=True,compact_nodes=True)
        
        
        for trial in range(start_trial,N_trials ):
            startime = timeit.default_timer()
            kd_tree_pathid = {}
            for count,path_id in enumerate(path_id_list_sample):
                print('Point_cloud, trial:',trial,'path count:',count,len_list)
                #startime = timeit.default_timer()
                segment_list, radius_list = randomize_segments_with_thickness(segments_path_dict[path_id],radius_path_dict[path_id],trial*(path_id+1) + count)
                thick_line_list = []
                for i in range(len(segment_list)):
                     if not(np.array_equal(segment_list[i][0],segment_list[i][1])):
                         thick_line_list += create_thick_line_point_list(segment_list[i][0],segment_list[i][1],segment_resolution,t*radius_list[i][0],t*radius_list[i][1])
                kd_tree_pathid[path_id] = KDTree(thick_line_list,balanced_tree=True,compact_nodes=True)
                #endtime = timeit.default_timer()
                #print('Time for path',count,'/',len(path_id_list) ,'is',endtime - startime)

            print('KD prepartion for trial',trial)
            kd_tree_radius = np.min([segment_distance_mean,t*radius_mean_1,t*radius_mean_2 ])
            

            print("Starting distance computation for trial",trial)
            intersections_paths = {}
            #kd_tree_radius = 0.005
            for count, path_id in enumerate(path_id_list_sample):
                #startime = timeit.default_timer()
                print('Distance computation, trial:',trial,'path count:',count,len_list)
                path_intersection_dict = {}

                intersections =  kd_tree_pathid[path_id].query_ball_tree(kd_tree_large ,r=kd_tree_radius)
                bodyid_intersected = intersecting_bodyid_from_kd_tree(intersections,dict_of_all_points)
                final_bodyid_list = []
                for bodyid in bodyid_intersected:
                        final_bodyid_list.append(bodyid)
                intersections_paths[path_id] =  list(set(final_bodyid_list +  adjacent_paths[path_id]))
                intersections_paths[path_id].remove(path_id)
                #intersections_paths[path_id] =  final_bodyid_list
            t_dict_results[t][trial] = intersections_paths
            endtime = timeit.default_timer()
            print('Parameter t',t,'time for trial',trial,endtime - startime)

            #with open(name +"_t_" +str(t) + "_path_intersections_kd_tree_radius_" +str(kd_tree_radius) + "_trial_" + str(trial) + ".pkl", "wb") as h:
            #    pickle.dump(intersections_paths , h)
            with open(path_save + name + '_directed_metagraph_dict_results.pkl', "wb") as h:
                    pickle.dump(t_dict_results, h)

