import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
import ast
import seaborn as sns
from scipy.stats import linregress
import math
from scipy.spatial import distance as dst
from scipy.spatial import ConvexHull

def create_line_point_list(start_point, end_point, density):

    parameter_list = np.linspace(0,1,density)
    line_point_list = []

    for t in parameter_list:
            x = (end_point[0] - start_point[0]) * t + start_point[0]
            y = (end_point[1] - start_point[1]) * t + start_point[1]
            z = (end_point[2] - start_point[2]) * t + start_point[2]
            line_point_list.append([x, y, z])
def create_line_point_radii_list(start_point, end_point,radius_start,radius_end,density):

    parameter_list = np.linspace(0,1,density)
    line_point_list = []

    for t in parameter_list:
            x = (end_point[0] - start_point[0]) * t + start_point[0]
            y = (end_point[1] - start_point[1]) * t + start_point[1]
            z = (end_point[2] - start_point[2]) * t + start_point[2]
            radius = (radius_end - radius_start) * t + radius_start
            line_point_list.append([x, y, z,radius])

    return line_point_list

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
        if 2 * radius_2 > distance_between_circles:
            number_of_inner_circles = np.max([math.ceil( 2 * radius_2 / distance_between_circles),3])
            inner_radius_list = [radius_2 * parameter for parameter in np.linspace(0,1, number_of_inner_circles)[1:-1]]
            for inner_radius in inner_radius_list:
                center_point = orientation_point
                radius = inner_radius
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

        return thick_points_list + line_point_list 
    
def bounding_box_from_merged_skeleton(merged_skeleton_dataset):
    x_min = np.min([np.min(merged_skeleton_dataset['x_1']),np.min(merged_skeleton_dataset['x_2'])])
    x_max = np.max([np.max(merged_skeleton_dataset['x_1']),np.max(merged_skeleton_dataset['x_2'])])

    y_min = np.min([np.min(merged_skeleton_dataset['y_1']),np.min(merged_skeleton_dataset['y_2'])])
    y_max = np.max([np.max(merged_skeleton_dataset['y_1']),np.max(merged_skeleton_dataset['y_2'])])

    z_min = np.min([np.min(merged_skeleton_dataset['z_1']),np.min(merged_skeleton_dataset['z_2'])])
    z_max = np.max([np.max(merged_skeleton_dataset['z_1']),np.max(merged_skeleton_dataset['z_2'])])

    return [x_min,x_max,y_min,y_max,z_min,z_max]

name_list = name_list = ['root_1','root_2','human_neuron','zebrafish_neuron','monkey_neuron',
                'rat_neuron','anthill','vascular_1','vascular_2','vascular_3',
                'mitochondrial','fruit_fly_1','fruit_fly_2','fruit_fly_3','fruit_fly_4']
final_results ={}
for name in name_list:
    print('**** Network:',name)
    print('Starting computation')
    final_results[name] = {}
    path = '../1. data/3. final_data/'
    skeleton_paths = pd.read_csv(path +  name+ '.paths.csv',index_col=[0])
    skeleton_paths.drop_duplicates(inplace=True)
    bounds = bounding_box_from_merged_skeleton(skeleton_paths)
    bounding_box_volume = abs(bounds[0] - bounds[1]) * abs(bounds[2] - bounds[3]) * abs(bounds[4] - bounds[5])
    longest_side = np.max([abs(bounds[0] - bounds[1]),abs(bounds[2] - bounds[3]),abs(bounds[4] - bounds[5])])
    shortest_side = np.min([abs(bounds[0] - bounds[1]),abs(bounds[2] - bounds[3]),abs(bounds[4] - bounds[5])])
    skeleton_paths['distance'] = np.sqrt((skeleton_paths['x_2']-  skeleton_paths['x_1'])**2 +  
                                               (skeleton_paths['y_2']-  skeleton_paths['y_1'])**2 +
                                               (skeleton_paths['z_2']-  skeleton_paths['z_1'])**2 )
    cube_side =  shortest_side *0.32
    cube_side_list = [proportion / 100 * cube_side for proportion in np.arange(15, 85, 5)]
    segment_resolution = np.min([skeleton_paths['distance'].mean(),skeleton_paths['radius_1'].mean(),skeleton_paths['radius_2'].mean()])/2
    final_results[name]['cube_side_list'] = cube_side_list
    
    print('Skeleton multiplication')
    multiplied_skeleton_part_merged = []
    for row in skeleton_paths[['pt_id_1','x_1','y_1','z_1','radius_1','pt_id_2','x_2','y_2','z_2','radius_2']].values:
        start_point = row[1:4]
        radius_start = row[4]
        end_point = row[6:9]
        radius_end = row[9]
        point_to_point_dist = dst.euclidean(start_point,end_point)
        ratio = (point_to_point_dist+radius_start+radius_end)/ (cube_side)
        #if ratio > 0.2:
        #    multiplication_factor =  7
        #    new_points = create_line_point_radii_list(start_point, end_point, radius_start,radius_end,multiplication_factor)
        #    for i in range(0,multiplication_factor-1):
        #        multiplied_skeleton_part_merged += [list(new_points[i])+list(new_points[i+1])]
        #else:
        new_row = [list(start_point) + [radius_start] + list(end_point) + [radius_end]]
        multiplied_skeleton_part_merged += new_row 

    multiplied_skeleton = pd.DataFrame(multiplied_skeleton_part_merged,columns=['x_1','y_1','z_1','radius_1','x_2','y_2','z_2','radius_2'])
    skeleton_paths = multiplied_skeleton.copy()
    print('Skeleton multiplication complete')
    
    point_cloud_list = []
    for i,row in enumerate(skeleton_paths[['x_1','y_1','z_1','radius_1','x_2','y_2','z_2','radius_2']].values):
        print('**** Network:',name)
        print( 'Computing point clouds',i,len(skeleton_paths))
        start_point = row[0:3]
        radius_start = row[3]
        end_point = row[4:7]
        radius_end = row[7]
        point_cloud_list += create_thick_line_point_list(start_point,end_point,segment_resolution, radius_start,radius_end)
    point_cloud_skeleton_paths  = pd.DataFrame(point_cloud_list,columns=['x','y','z'])
    filled_boxes_count_list = []
    for cube_side in cube_side_list:
        print('**** Network:',name,'Cube side',cube_side,'out of',cube_side_list)
        print('Cube side ratio:',cube_side/longest_side)
        x_axis = np.arange(bounds[0],bounds[1]+cube_side,cube_side)
        y_axis = np.arange(bounds[2],bounds[3]+cube_side,cube_side)
        z_axis = np.arange(bounds[4],bounds[5]+cube_side,cube_side)
        number_of_segments_used = 0
        filled_boxes_count = 0
        empty_boxes_count = 0
        box_count = 0
        cube_dict = {}
        for i in range(len(x_axis)-1):
                print(i,len(x_axis))
                box_count += 1
                for j in range(len(y_axis)-1):
                            for k in range(len(z_axis)-1):
                                points_inside = point_cloud_skeleton_paths[(point_cloud_skeleton_paths['x'] >= x_axis[i]) & 
                                                            (point_cloud_skeleton_paths['x'] <= x_axis[i+1]) &
                                                            (point_cloud_skeleton_paths['y'] >= y_axis[j]) &
                                                            (point_cloud_skeleton_paths['y'] <= y_axis[j+1]) &
                                                            (point_cloud_skeleton_paths['z'] >= z_axis[k]) &
                                                            (point_cloud_skeleton_paths['z'] <= z_axis[k+1])]
                                if len(points_inside) == 0:
                                    empty_boxes_count += 1
                                else:
                                    filled_boxes_count += 1
        filled_boxes_count_list.append(filled_boxes_count)
    final_results[name]['filled_boxes_count_list'] = filled_boxes_count_list
    with open("2. fractal_dimension_results/"+name + "_basic_results.pkl", "wb") as h:
           pickle.dump(final_results[name], h)    
