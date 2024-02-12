import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import random
import pandas as pd
import pickle

plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rc('xtick', labelsize=16) 
plt.rc('ytick', labelsize=16) 
def e_distance(start_point,end_point):
    distance = np.sqrt((end_point[0] - start_point[0])**2 + (end_point[1] - start_point[1] )**2 + (end_point[2] - start_point[2] )**2)
    return distance

def rw_length(rw):
    total_distance = 0
    for i in range(len(rw)-1):
        start = rw[i][0:3]
        end = rw[i+1][0:3]
        total_distance += distance.euclidean(start,end)
        
    return np.round(total_distance,4)

name_list = ['fruit_fly_1','fruit_fly_2','fruit_fly_3','fruit_fly_4']
clustering_steps_dict = {}
for name in name_list:
    print('*****',name)
    path_source = '5. intersection_connectome/'
    path_save = '6. intersection_connectome_skeleton_resolution/'
    df = pd.read_csv(path_source + name + '.paths.csv',index_col = [0])
    #stopping_criteria_search = True
    new_path_dict = {}
    stopping_step_dict = {}
    path_id_list = list(set(df['path_id'].values))
    for count,path_id in enumerate(path_id_list):
        print(count,len(path_id_list))
        df_one_path = df[df['path_id']==path_id]
        bodyid_pre = df_one_path['source'].values[0]
        bodyid_post = df_one_path['target'].values[0]
        rw = []
        #rw.append(df_one_path[['x_1','y_1','z_1','pt_id_1','radius_1']].values[0])
        for point in df_one_path[['x_1','y_1','z_1','pt_id_1','radius_1']].values:
            rw.append(np.array(point))
        rw.append(df_one_path[['x_2','y_2','z_2','pt_id_2','radius_2']].values[-1])

        total_path_length = rw_length(rw)
        mean_cos_angle_list =  []
        std_cos_angle_list = []
        removed_cos_angle_list = []
        #print('Number of segments after correction',len(rw))
        segment_list = []

        length_ratio_list = []
        #Looks for duplicates or nan vectors (length 0)
        for index,segment in enumerate(rw):
            if index <len(rw)-1:
                start = list(rw[index][0:3])
                end = list(rw[index+1][0:3])
                diff_vector = np.array(end) - np.array(start)
                start_part = list(rw[index])
                end_part = list(rw[index+1])
                if (np.linalg.norm(np.array(diff_vector))!= 0):
                        diff_vector_hat = np.array(diff_vector) / np.linalg.norm(np.array(diff_vector))
                        segment_list.append(start_part + end_part + list(diff_vector_hat))

        clustering_steps = len(segment_list)-2
        clustering_steps_list = list(range(0,clustering_steps+1))  
        stopping_criteria_search = True
        for step in range(0,clustering_steps+2):
            print('Clustering step',step,'out of',clustering_steps+1)
            segment_list = []
            for index,segment in enumerate(rw):
                if index <len(rw)-1:
                    start = list(rw[index][0:3])
                    end = list(rw[index+1][0:3])
                    diff_vector = np.array(end) - np.array(start)
                    diff_vector_hat = np.array(diff_vector) / np.linalg.norm(np.array(diff_vector))
                    segment_list.append(list(rw[index]) +list(rw[index+1]) + list(diff_vector_hat))
                    #if step == 0:
                    #    if (np.linalg.norm(np.array(diff_vector))!= 0):
                    #        diff_vector_hat = np.array(diff_vector) / np.linalg.norm(np.array(diff_vector))
                    #        segment_list.append(list(rw[index]) + list(rw[index+1]) + list(diff_vector_hat))
                    #else:
                    #    diff_vector_hat = np.array(diff_vector) / np.linalg.norm(np.array(diff_vector))
                    #    segment_list.append(list(rw[index]) +list(rw[index+1]) + list(diff_vector_hat))
            df_segments = pd.DataFrame(segment_list,columns = ['x_1','y_1','z_1','pt_id_1','radius_1','x_2','y_2','z_2','pt_id_2','radius_2','x_norm','y_norm','z_norm'])
            direction_vector_list = df_segments[['x_norm','y_norm','z_norm']].values
            cos_angle_dict = {}
            for i,norm_diff_vec in enumerate(direction_vector_list):
                if i <len(direction_vector_list)-1:
                    cos_angle = np.dot(direction_vector_list[i],direction_vector_list[i+1])
                    if i == 0 or i == len(direction_vector_list)-2:
                     # Fixing the starting and ending segment, as to avoiding changing initial node positions.
                        cos_angle_dict[(i,i+1)] = 0
                    else:
                        cos_angle_dict[(i,i+1)] = cos_angle
            sorted_segment_indices = sorted(cos_angle_dict, key=cos_angle_dict.get, reverse=True)
            #print(sorted_segment_indices)
            if len(cos_angle_dict) > 0:
                mean_cos_angle_list.append(np.mean(list(cos_angle_dict.values())))
                std_cos_angle_list.append(np.std(list(cos_angle_dict.values())))
                removed_cos_angle_list.append(cos_angle_dict[sorted_segment_indices[0]] )
                #print('mean',np.mean(list(cos_angle_dict.values())),'std',np.std(list(cos_angle_dict.values())))
                #Dropping
                #print(df_segments)
                index_1 = sorted_segment_indices[0][0]
                index_2 =  sorted_segment_indices[0][1]
                #print(index_1,index_2)
                first_row = df_segments.iloc[index_1 ]
                second_row = df_segments.iloc[ index_2]
                new_start = first_row[['x_1','y_1','z_1','pt_id_1','radius_1']].values
                new_end = second_row[['x_2','y_2','z_2','pt_id_2','radius_2']].values
                new_diff_vec =np.array(new_end[0:3]) - np.array(new_start[0:3])
                new_diff_vec_hat = np.array(diff_vector) / np.linalg.norm(np.array(diff_vector))
                df_dropped = df_segments.drop(index=[index_1 , index_2])
                #Inserting
                df_dropped.loc[(index_1 + index_2)/2] = new_start[0], new_start[1],new_start[2], new_start[3], new_start[4],new_end[0],new_end[1],new_end[2],new_end[3],new_end[4],new_diff_vec_hat[0],new_diff_vec_hat[1],new_diff_vec_hat[2]
            else:
                df_dropped = df_segments.copy()
            df_reset = df_dropped.sort_index().reset_index(drop=True)
            rw_new = []
            for point in df_reset[['x_1','y_1','z_1','pt_id_1','radius_1']].values:
                rw_new.append(np.array(point))
            rw_new.append(df_reset[['x_2','y_2','z_2','pt_id_2','radius_2']].values[-1])
            length_ratio = rw_length(rw_new)/total_path_length 
            length_ratio_list.append(length_ratio)
            if length_ratio< 0.9999999 and stopping_criteria_search == True:
                stopping_step = step
                stopping_criteria_search = False
                stopping_step_dict[path_id] = step
                break
            rw = rw_new
        if stopping_criteria_search == True:
            stopping_step_dict[path_id] = step
        new_path_dict[(path_id,bodyid_pre,bodyid_post)] = rw
    clustering_steps_dict[name] = stopping_step_dict
    row_list = []
    for path_pre_post in new_path_dict.keys():
        point_list = new_path_dict[path_pre_post]
        for i in range(len(point_list)-1):
            #print(i)
            point_1 = point_list[i]
            point_2 = point_list[i+1]
            coord_1 =  point_list[i][0:3]
            coord_2 =  point_list[i+1][0:3]
            pt_id_1 = point_list[i][3]
            pt_id_2 =point_list[i+1][3]
            if type(pt_id_1) != str:
                pt_id_1 = str(int(pt_id_1))
            if type(pt_id_2) != str:
                pt_id_2 = str(int(pt_id_2))
            radius_1 =point_list[i][4]
            radius_2 =point_list[i+1][4]
            row = list( coord_1) + [pt_id_1] + [radius_1] +list( coord_2) + [pt_id_2] + [radius_2] +[path_pre_post[0]] + [path_pre_post[1]] + [path_pre_post[2]]
            row_list.append(row)
    new_df = pd.DataFrame(row_list,columns=['x_1','y_1','z_1','pt_id_1','radius_1','x_2','y_2','z_2','pt_id_2','radius_2','path_id','source','target'])
    new_df = new_df[['pt_id_1','x_1','y_1','z_1','radius_1','pt_id_2','x_2','y_2','z_2','radius_2','path_id','source','target']].copy()
    new_df.to_csv(path_save + name + '.paths.csv')
    

print("Checking if paths start and end at bodyids")
for name in name_list:
    print('*****',name)
    path_save = '6. intersection_connectome_skeleton_resolution/'
    df = pd.read_csv(path_save + name + '.paths.csv',index_col = [0])
    for path_id in set(list(df['path_id'].values)):
        bool_1 = df[df['path_id'] == path_id]['pt_id_1'].values[0] == df[df['path_id'] == path_id]['source'].values[0]
        bool_2 = df[df['path_id'] == path_id]['pt_id_2'].values[-1]== df[df['path_id'] == path_id]['target'].values[0]
        if bool_1 == False:
            print('False! Clustering messed up the connectome pt_id_1 for ',name)
            break
        if bool_2 == False:
            print('False! Clustering messed up the connectome pt_id_2 for ',name)
            break

print('Checking minimum distance')
for name in name_list:
    print('*****',name)
    path_save = '6. intersection_connectome_skeleton_resolution/'
    df = pd.read_csv(path_save + name + '.paths.csv',index_col = [0])
    df['distance'] = ((df['x_2']-df['x_1'])**2 + (df['y_2']-df['y_1'])**2  + (df['z_2']-df['z_1'])**2)**1/2
    df = df[df['distance'] != 0].copy()
    print("The shortest segment is of length",np.min(df['distance']))
    df.to_csv(path_save + name + '.paths.csv')