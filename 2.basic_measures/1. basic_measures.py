import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
import ast
import seaborn as sns
from scipy.stats import linregress
#sns.set_theme(style="whitegrid")
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

#name_list = [ 'vascular_2','vascular_3','vascular_1','human_neuron','rat_neuron','zebrafish_neuron','mitochondrial','anthill','root_1','root_2','monkey_neuron','fruit_fly_2','fruit_fly_3','fruit_fly_1','fruit_fly_4']
#name_list = ['human_neuron','rat_neuron',
#            'monkey_neuron','zebrafish_neuron',
#    'vascular_2','vascular_3','vascular_1','mitochondrial','anthill','root_1','root_2',
name_list = ['fruit_fly_1','fruit_fly_4']
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
    skeleton_paths['volume'] = 1/3 * np.pi * (skeleton_paths['radius_1']**2 + skeleton_paths['radius_2']**2 +    skeleton_paths['radius_1']*skeleton_paths['radius_2'])*skeleton_paths['distance']
    skeleton_paths_original = skeleton_paths.copy()
    final_results[name]['mean_segment_length'] = np.mean(skeleton_paths['distance'])
    final_results[name]['segment_length_list'] = skeleton_paths['distance'].values
    final_results[name]['mean_radius_length'] = np.mean(skeleton_paths['radius_1'])
    final_results[name]['radius_list'] = skeleton_paths['radius_1'].values
    final_results[name]['mean_segment_over_bounding_box_average_length'] = np.mean(skeleton_paths['distance'])/(np.mean([abs(bounds[0] - bounds[1]),
                                                                                                                          abs(bounds[2] - bounds[3]),
                                                                                                                          abs(bounds[4] - bounds[5])]))
    cube_side =  shortest_side *0.32
    cube_side_list = [proportion/100 * cube_side for proportion in np.arange(15,105,5)]
    segment_resolution = np.min([skeleton_paths['distance'].mean(),skeleton_paths['radius_1'].mean(),skeleton_paths['radius_2'].mean()])/2
    final_results[name]['cube_side'] = cube_side
    final_results[name]['cube_side_list'] = cube_side_list
    final_results[name]['segment_resolution_distance'] = segment_resolution
    
    print('Skeleton multiplication')
    multiplied_skeleton_part_merged = []
    for row in skeleton_paths[['pt_id_1','x_1','y_1','z_1','radius_1','pt_id_2','x_2','y_2','z_2','radius_2']].values:
        start_point = row[1:4]
        radius_start = row[4]
        end_point = row[6:9]
        radius_end = row[9]
        point_to_point_dist = dst.euclidean(start_point,end_point)
        ratio = (point_to_point_dist+radius_start+radius_end)/ (cube_side)
        if ratio > 0.2:
            multiplication_factor = 2 + 5
           # print('Multiplying skeleton',multiplying_factor )
            new_points = create_line_point_radii_list(start_point, end_point, radius_start,radius_end,multiplication_factor)
            for i in range(0,multiplication_factor-1):
                multiplied_skeleton_part_merged += [list(new_points[i])+list(new_points[i+1])]
        else:
            new_row = [list(start_point) + [radius_start] + list(end_point) + [radius_end]]
            multiplied_skeleton_part_merged += new_row 

    multiplied_skeleton = pd.DataFrame(multiplied_skeleton_part_merged,columns=['x_1','y_1','z_1','radius_1','x_2','y_2','z_2','radius_2'])
    skeleton_paths = multiplied_skeleton.copy()
    print('Skeleton multiplication complete')
    
    point_cloud_list = []
    for i,row in enumerate(skeleton_paths[['x_1','y_1','z_1','radius_1','x_2','y_2','z_2','radius_2']].values):
        print('**** Network:',name)
        print( 'Computing fractal dimension',i,len(skeleton_paths))
        start_point = row[0:3]
        radius_start = row[3]
        end_point = row[4:7]
        radius_end = row[7]
        point_cloud_list += create_thick_line_point_list(start_point,end_point,segment_resolution, radius_start,radius_end)
    point_cloud_skeleton_paths  = pd.DataFrame(point_cloud_list,columns=['x','y','z'])
    link_paths = pd.read_csv(path +  name+ '.paths.csv',index_col=[0])
    link_paths = link_paths[['x_1','y_1','z_1','x_2','y_2','z_2','path_id']].copy()
    filled_boxes_count_list = []
    all_densities = []
    all_link_densities = []
    number_of_segments_used_list = []
    for cube_side in cube_side_list:
        print('**** Network:',name,'Cube side',cube_side,'out of',cube_side_list)
        print('Cube side ratio:',cube_side/longest_side)
        x_axis = np.arange(bounds[0],bounds[1],cube_side)
        y_axis = np.arange(bounds[2],bounds[3],cube_side)
        z_axis = np.arange(bounds[4],bounds[5],cube_side)
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
                                points_inside = point_cloud_skeleton_paths[(point_cloud_skeleton_paths['x'] > x_axis[i]) & 
                                                            (point_cloud_skeleton_paths['x'] < x_axis[i+1]) &
                                                            (point_cloud_skeleton_paths['y'] > y_axis[j]) &
                                                            (point_cloud_skeleton_paths['y'] < y_axis[j+1]) &
                                                            (point_cloud_skeleton_paths['z'] > z_axis[k]) &
                                                            (point_cloud_skeleton_paths['z'] < z_axis[k+1])]
                                links_inside  = link_paths[(link_paths['x_1'] > x_axis[i]) & 
                                                            (link_paths['x_1'] < x_axis[i+1]) &
                                                            (link_paths['y_1'] > y_axis[j]) &
                                                            (link_paths['y_1'] < y_axis[j+1]) &
                                                            (link_paths['z_1'] > z_axis[k]) &
                                                            (link_paths['z_2'] < z_axis[k+1])]
                                if len(links_inside) == 0:
                                    all_link_densities.append(0)
                                else:
                                    all_link_densities.append(len(set(links_inside['path_id'].values)))
                                if len(points_inside) == 0:
                                    empty_boxes_count += 1
                                else:
                                    filled_boxes_count += 1
                                    if cube_side == cube_side_list[0] and len(points_inside)>=4:
                                        convex_hull = ConvexHull(points_inside.values,qhull_options="QJ")
                                        local_density = (convex_hull.volume)/cube_side**3
                                        all_densities.append(local_density)
                                        if local_density > 1 :
                                            print("Density larger than 1")
                                            break
        filled_boxes_count_list.append(filled_boxes_count)
    fractal_dimension =  -linregress(np.log(cube_side_list),np.log(filled_boxes_count_list))[0]
    final_results[name]['filled_boxes_count_list'] = filled_boxes_count_list
    final_results[name]['fractal_dimension'] = np.round(fractal_dimension,3)
    final_results[name]['density_list'] = all_densities
    final_results[name]['density_max'] = np.max(all_densities)
    final_results[name]['link_density_list'] =all_link_densities
    with open("1. results_basic_properties/"+name + "_basic_results.pkl", "wb") as h:
           pickle.dump(final_results[name], h)    

            
    infile = open('1. results_basic_properties/' + name +'_basic_results.pkl','rb')
    final_results[name] = pickle.load(infile)

    print('Computing complementary straightness')
    skeleton_paths = pd.read_csv(path +  name+ '.paths.csv',index_col=[0])
    skeleton_paths.drop_duplicates(inplace=True)
    skeleton_paths['distance'] = np.sqrt((skeleton_paths['x_2']-  skeleton_paths['x_1'])**2 +  
                                               (skeleton_paths['y_2']-  skeleton_paths['y_1'])**2 +
                                               (skeleton_paths['z_2']-  skeleton_paths['z_1'])**2 )
    skeleton_paths['volume'] = 1/3 * np.pi * (skeleton_paths['radius_1']**2 + skeleton_paths['radius_2']**2 +    skeleton_paths['radius_1']*skeleton_paths['radius_2'])*skeleton_paths['distance']
    segments_path_id = skeleton_paths[['bodyId_pre','bodyId_post','path_id']].copy()
    segments_path_id.drop_duplicates(inplace=True) 
    
    connectome = segments_path_id
    connectome_list = connectome.values.tolist()
    connectome_loc = []
    for i,pair in enumerate(connectome_list):
        #print(i,pair)
        start_row = skeleton_paths[(skeleton_paths['pt_id_1'] == pair[0] ) & (skeleton_paths['path_id'] == pair[2] )]
        if len(start_row) == 0:
            start_row = skeleton_paths[(skeleton_paths['pt_id_2'] == pair[0] ) & (skeleton_paths['path_id'] == pair[2] )]
            start_loc = [start_row['x_2'].values[0],start_row['y_2'].values[0],start_row['z_2'].values[0]]
        else:
            start_loc = [start_row['x_1'].values[0],start_row['y_1'].values[0],start_row['z_1'].values[0]]

        end_row = skeleton_paths[(skeleton_paths['pt_id_2'] == pair[1] )& (skeleton_paths['path_id'] == pair[2] )]
        if len(end_row) == 0:
            end_row = skeleton_paths[(skeleton_paths['pt_id_1'] == pair[1] )& (skeleton_paths['path_id'] == pair[2] )]
            end_loc = [end_row['x_1'].values[0],end_row['y_1'].values[0],end_row['z_1'].values[0]]
        else:
            end_loc = [end_row['x_2'].values[0],end_row['y_2'].values[0],end_row['z_2'].values[0]]
        pair.append(start_loc)
        pair.append(end_loc)
        connectome_loc.append(pair)
        
    euclidean_path_length_list = []
    total_path_length_list = []
    contraction_list = []
    count = 0
    contraction_pairs = {}
    volume_pairs = {}
    link_path_length_dict = {}
    nr_segments_per_link_dict = {}
    for index,pair in enumerate(connectome_loc):
        bodyid_pre = pair[0]
        bodyid_post = pair[1]
        skeleton_df = skeleton_paths[skeleton_paths['path_id'] == pair[2]].copy()
        nr_segments_per_link = len(skeleton_df)
        point_list = skeleton_df[['x_1','y_1','z_1']].values.tolist()
        radius_list = skeleton_df[['radius_1']].values.tolist()
        distance_list = []
        volume_list =  []
        for i,row in skeleton_df.iterrows():
            point_start =  row[['x_1','y_1','z_1']].values.tolist()
            point_end =  row[['x_2','y_2','z_2']].values.tolist()
            radius_start =  row[['radius_1']].values.tolist()[0]
            radius_end = row[['radius_2']].values.tolist()[0]
            distance = np.sqrt((point_end[0] - point_start[0])**2 + (point_end[1] - point_start[1])**2 + (point_end[2] - point_start[2])**2 )
            distance_list.append(distance)
            volume = 1/3 * np.pi * distance * ( radius_start**2 + radius_end**2 + radius_start*radius_end)
            volume_list.append(volume)
        total_path_length = np.sum(distance_list)
        total_volume = np.sum(volume_list)
        if total_path_length > 0:
            total_path_length_list.append(total_path_length)
            node_source = pair[3]
            node_target = pair[4]
            euclidean_distance = np.sqrt((node_source[0]-node_target[0])**2 +(node_source[1]-node_target[1])**2 + (node_source[2]-node_target[2])**2) 
            euclidean_path_length_list.append(euclidean_distance)
            contraction = euclidean_distance/total_path_length
            contraction_list.append(contraction )
            contraction_pairs[pair[2]] = 1 - contraction
            volume_pairs[pair[2]] = total_volume
            link_path_length_dict[pair[2]] = total_path_length
            nr_segments_per_link_dict[pair[2]] = nr_segments_per_link
            if contraction > 1:
                print('This value cannot be larger than 1')
                break
        count +=1 
    total_volume_of_the_network = np.sum(list(volume_pairs.values()))
    final_results[name]['link_volume_normed'] = {}
    for link_volume_pair in volume_pairs.keys():
        final_results[name]['link_volume_normed'][link_volume_pair] = final_results[name]['link_volume'][link_volume_pair]/total_volume_of_the_network
    final_results[name]['c_straightness'] = contraction_pairs
    final_results[name]['link_volume'] = volume_pairs
    final_results[name]['link_segments'] = nr_segments_per_link_dict
    final_results[name]['link_path_length'] = link_path_length_dict
    with open("1. results_basic_properties_new/"+name + "_basic_results.pkl", "wb") as h:
               pickle.dump(final_results[name], h)