# -*- coding: utf-8 -*-
"""
Operator Index for metamaterial design problems with Bootstrapping

@author: roshan94
"""
from pygmo import hypervolume
import csv
#import math
import numpy as np
from IPython.core.debugger import set_trace

### Useful functions
def get_fileloc(truss_prob, model, run_mode, run_num):
    filepath_init = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\Operator Index Data\\' # for office system
    #filepath_init = 'C:\\Users\\rosha\\Documents\\SEAK Lab Github\\KD3M3\\result\\Operator Index Data\\' # for home system

    if run_mode == 1: 
        filepath_mode = 'Random\\'
        filename_init = 'random_operator_index_data_'
    elif run_mode == 2:
        filename_init = 'EpsilonMOEA_emoea__'
        filepath_mode = 'Epsilon MOEA\\'

    filepath_problem = 'Artery Problem\\'
    filename_problem = 'artery_'
    
    if truss_prob:
        filepath_problem = 'Truss Problem\\'
        filename_problem = 'prob2_'
    
    if model == 1:
        filepath_model = 'Fiber Model\\'
        filename_model = 'fibre'
    elif model == 2:
        filepath_model = 'Truss Model\\'
        filename_model = 'truss'
    elif model == 3:
        filepath_model = 'Beam Model\\'
        filename_model = 'beam'
    
    filename = filename_init + filename_problem + filename_model + str(run_num) + '.csv'
    filepath = filepath_init + filepath_problem + filepath_model + filepath_mode + filename
    
    return filepath

def compute_pareto_front(population):
    pop_size = len(population)
    obj_num = 2

    domination_counter = [0] * pop_size

    for i in range(pop_size):
        for j in range(i+1, pop_size):
            # check each objective for dominance
            dominate = [0] * obj_num
            for k in range(obj_num):
                if population[i][k] > population[j][k]:
                    dominate[k] = 1
                elif population[i][k] < population[j][k]:
                    dominate[k] = -1
            if -1 not in dominate and 1 in dominate:
                domination_counter[i] += 1
            elif -1 in dominate and 1 not in dominate:
                domination_counter[j] += 1

    pareto_solutions = []
    for i in range(len(domination_counter)):
        if domination_counter[i] == 0:
            pareto_solutions.append(population[i])
    return pareto_solutions

def get_pen_objs(norm_objs_array, feas_array, conn_array, stiffrat_array):
    feas_viol_array = np.subtract(1, feas_array)
    conn_viol_array = np.subtract(1, conn_array)
    pen_weight = np.ones((3))
    
    feas_pen_add = np.multiply(feas_viol_array, pen_weight[0])
    conn_pen_add = np.multiply(conn_viol_array, pen_weight[1])
    stiffrat_pen_add = np.multiply(stiffrat_array, pen_weight[2])
    pen_total_add1 = np.add(feas_pen_add, conn_pen_add)
    pen_total_add = np.add(pen_total_add1, stiffrat_pen_add)
    
    norm_objs1_array = [x[0] for x in norm_objs_array]
    norm_objs2_array = [x[1] for x in norm_objs_array]
    
    pen_objs1_array = np.add(norm_objs1_array, pen_total_add)
    pen_objs2_array = np.add(norm_objs2_array, pen_total_add)
    
    return np.column_stack((pen_objs1_array, pen_objs2_array))
    

def compute_hv(population):
    array_archs = np.zeros((len(population), 2))
    for i in range(len(population)):
        array_archs[i] = population[i]
    hv_object = hypervolume(array_archs)
    hv = hv_object.compute([1.1,1.1])/1.1**2
    return hv

def read_csv(fileloc_csv, truss_prob):
    
    with open(fileloc_csv, newline='') as csvfile:
        data = [row for row in csv.reader(csvfile)]
        
        norm_objs1_dat = np.zeros(len(data)-1)
        norm_objs2_dat = np.zeros(len(data)-1)
        feas_dat = np.zeros(len(data)-1)
        conn_dat = np.zeros(len(data)-1)
        
        norm_objs1_partcoll_dat = np.zeros(len(data)-1)
        norm_objs2_partcoll_dat = np.zeros(len(data)-1)
        feas_partcoll_dat = np.zeros(len(data)-1)
        conn_partcoll_dat = np.zeros(len(data)-1)
        
        norm_objs1_nodalprop_dat = np.zeros(len(data)-1)
        norm_objs2_nodalprop_dat = np.zeros(len(data)-1)
        feas_nodalprop_dat = np.zeros(len(data)-1)
        conn_nodalprop_dat = np.zeros(len(data)-1)
        
        norm_objs1_orient_dat = np.zeros(len(data)-1)
        norm_objs2_orient_dat = np.zeros(len(data)-1)
        feas_orient_dat = np.zeros(len(data)-1)
        conn_orient_dat = np.zeros(len(data)-1)
        
        norm_objs1_inters_dat = np.zeros(len(data)-1)
        norm_objs2_inters_dat = np.zeros(len(data)-1)
        feas_inters_dat = np.zeros(len(data)-1)
        conn_inters_dat = np.zeros(len(data)-1)
        
        if truss_prob:
            stiffrat_dat = np.zeros(len(data)-1)
            stiffrat_partcoll_dat = np.zeros(len(data)-1)
            stiffrat_nodalprop_dat = np.zeros(len(data)-1)
            stiffrat_orient_dat = np.zeros(len(data)-1)
            stiffrat_inters_dat = np.zeros(len(data)-1)
        
        valid_count = 0
        for x in range(len(data)-1):
            if truss_prob:
                data_float = list(map(float,data[x+1][1:6]))
                data_partcoll_float = list(map(float,data[x+1][7:12]))
                data_nodalprop_float = list(map(float,data[x+1][13:18]))
                data_orient_float = list(map(float,data[x+1][19:24]))
                data_inters_float = list(map(float,data[x+1][25:30]))
                
            else:
                data_float = list(map(float,data[x+1][1:5]))
                data_partcoll_float = list(map(float,data[x+1][6:10]))
                data_nodalprop_float = list(map(float,data[x+1][11:15]))
                data_orient_float = list(map(float,data[x+1][16:20]))
                data_inters_float = list(map(float,data[x+1][21:25]))
                
            if (any(np.isnan(np.array(data_float))) or any(np.isinf(np.array(data_float)))):
                continue
            
            if (any(np.isnan(np.array(data_partcoll_float))) or any(np.isinf(np.array(data_partcoll_float)))):
                continue
            
            if (any(np.isnan(np.array(data_nodalprop_float))) or any(np.isinf(np.array(data_nodalprop_float)))):
                continue
            
            if (any(np.isnan(np.array(data_orient_float))) or any(np.isinf(np.array(data_orient_float)))):
                continue
            
            if (any(np.isnan(np.array(data_inters_float))) or any(np.isinf(np.array(data_inters_float)))):
                continue
            
            norm_objs1_dat[valid_count] = data_float[0]
            norm_objs2_dat[valid_count] = data_float[1]
            feas_dat[valid_count] = data_float[2] 
            conn_dat[valid_count] = data_float[3]
            
            norm_objs1_partcoll_dat[valid_count] = data_partcoll_float[0]
            norm_objs2_partcoll_dat[valid_count] = data_partcoll_float[1]
            feas_partcoll_dat[valid_count] = data_partcoll_float[2]
            conn_partcoll_dat[valid_count] = data_partcoll_float[3]
            
            norm_objs1_nodalprop_dat[valid_count] = data_nodalprop_float[0]
            norm_objs2_nodalprop_dat[valid_count] = data_nodalprop_float[1]
            feas_nodalprop_dat[valid_count] = data_nodalprop_float[2] 
            conn_nodalprop_dat[valid_count] = data_nodalprop_float[3]
            
            norm_objs1_orient_dat[valid_count] = data_orient_float[0]
            norm_objs2_orient_dat[valid_count] = data_orient_float[1]
            feas_orient_dat[valid_count] = data_orient_float[2] 
            conn_orient_dat[valid_count] = data_orient_float[3]
            
            norm_objs1_inters_dat[valid_count] = data_inters_float[0]
            norm_objs2_inters_dat[valid_count] = data_inters_float[1]
            feas_inters_dat[valid_count] = data_inters_float[2] 
            conn_inters_dat[valid_count] = data_inters_float[3]
            
            if truss_prob:
                if (data_float[4] > 50):
                    stiffrat_dat[valid_count] = 5
                else:
                    stiffrat_dat[valid_count] = data_float[4]
                if (data_partcoll_float[4] > 50):
                    stiffrat_partcoll_dat[valid_count] = 5
                else:
                    stiffrat_partcoll_dat[valid_count] = data_partcoll_float[4]
                if (data_nodalprop_float[4] > 50):
                    stiffrat_nodalprop_dat[valid_count] = 5
                else:
                    stiffrat_nodalprop_dat[valid_count] = data_nodalprop_float[4]
                if (data_orient_float[4] > 50):
                    stiffrat_orient_dat[valid_count] = 5
                else:
                    stiffrat_orient_dat[valid_count] = data_orient_float[4]
                if (data_inters_float[4] > 50):
                    stiffrat_inters_dat[valid_count] = 5
                else:
                    stiffrat_inters_dat[valid_count] = data_inters_float[4]
                
            valid_count += 1
            
    norm_objs1 = norm_objs1_dat[:valid_count]
    norm_objs2 = norm_objs2_dat[:valid_count]
    feas = feas_dat[:valid_count]
    conn = conn_dat[:valid_count]
            
    norm_objs1_partcoll = norm_objs1_partcoll_dat[:valid_count]
    norm_objs2_partcoll = norm_objs2_partcoll_dat[:valid_count]
    feas_partcoll = feas_partcoll_dat[:valid_count]
    conn_partcoll = conn_partcoll_dat[:valid_count]
            
    norm_objs1_nodalprop = norm_objs1_nodalprop_dat[:valid_count]
    norm_objs2_nodalprop = norm_objs2_nodalprop_dat[:valid_count]
    feas_nodalprop = feas_nodalprop_dat[:valid_count]
    conn_nodalprop = conn_nodalprop_dat[:valid_count]
            
    norm_objs1_orient = norm_objs1_orient_dat[:valid_count]
    norm_objs2_orient = norm_objs2_orient_dat[:valid_count]
    feas_orient = feas_orient_dat[:valid_count]
    conn_orient = conn_orient_dat[:valid_count]
            
    norm_objs1_inters = norm_objs1_inters_dat[:valid_count]
    norm_objs2_inters = norm_objs2_inters_dat[:valid_count]
    feas_inters = feas_inters_dat[:valid_count]
    conn_inters = conn_inters_dat[:valid_count]
        
    if truss_prob:
        stiffrat = stiffrat_dat[:valid_count]
        stiffrat_partcoll = stiffrat_partcoll_dat[:valid_count]
        stiffrat_nodalprop = stiffrat_nodalprop_dat[:valid_count]
        stiffrat_orient = stiffrat_orient_dat[:valid_count]
        stiffrat_inters = stiffrat_inters_dat[:valid_count]
        
    norm_objs = np.column_stack((norm_objs1, norm_objs2))
    norm_objs_partcoll = np.column_stack((norm_objs1_partcoll, norm_objs2_partcoll))
    norm_objs_nodalprop = np.column_stack((norm_objs1_nodalprop, norm_objs2_nodalprop))
    norm_objs_orient = np.column_stack((norm_objs1_orient, norm_objs2_orient))
    norm_objs_inters = np.column_stack((norm_objs1_inters, norm_objs2_inters))
    
    data = {}
    
    data['objs'] = norm_objs
    data['objs_partcoll'] = norm_objs_partcoll
    data['objs_nodalprop'] = norm_objs_nodalprop
    data['objs_orient'] = norm_objs_orient
    data['objs_inters'] = norm_objs_inters
    
    data['feas'] = feas
    data['feas_partcoll'] = feas_partcoll
    data['feas_nodalprop'] = feas_nodalprop
    data['feas_orient'] = feas_orient
    data['feas_inters'] = feas_inters
    
    data['conn'] = conn
    data['conn_partcoll'] = conn_partcoll
    data['conn_nodalprop'] = conn_nodalprop
    data['conn_orient'] = conn_orient
    data['conn_inters'] = conn_inters
    
    if truss_prob:
        data['stiffrat'] = stiffrat
        data['stiffrat_partcoll'] = stiffrat_partcoll
        data['stiffrat_nodalprop'] = stiffrat_nodalprop
        data['stiffrat_orient'] = stiffrat_orient
        data['stiffrat_inters'] = stiffrat_inters
        
    return data

def get_pfs_dataset(data_rand, truss_prob):
    
    norm_objs_combined = data_rand['objs'].tolist()
    norm_objs_combined_partcoll = data_rand['objs_partcoll'].tolist()
    norm_objs_combined_nodalprop = data_rand['objs_nodalprop'].tolist()
    norm_objs_combined_orient = data_rand['objs_orient'].tolist()
    norm_objs_combined_inters = data_rand['objs_inters'].tolist()
    
    feas_combined = data_rand['feas'].tolist()
    feas_combined_partcoll = data_rand['feas_partcoll'].tolist()
    feas_combined_nodalprop = data_rand['feas_nodalprop'].tolist()
    feas_combined_orient = data_rand['feas_orient'].tolist()
    feas_combined_inters = data_rand['feas_inters'].tolist()
    
    conn_combined = data_rand['conn'].tolist()
    conn_combined_partcoll = data_rand['conn_partcoll'].tolist()
    conn_combined_nodalprop = data_rand['conn_nodalprop'].tolist()
    conn_combined_orient = data_rand['conn_orient'].tolist()
    conn_combined_inters = data_rand['conn_inters'].tolist()
    
    if truss_prob:
        stiffrat_combined = data_rand['stiffrat'].tolist()
        stiffrat_combined_partcoll = data_rand['stiffrat_partcoll'].tolist()
        stiffrat_combined_nodalprop = data_rand['stiffrat_nodalprop'].tolist()
        stiffrat_combined_orient = data_rand['stiffrat_orient'].tolist()
        stiffrat_combined_inters = data_rand['stiffrat_inters'].tolist()
    else:
        stiffrat_combined = np.zeros(len(feas_combined))
        stiffrat_combined_partcoll = np.zeros(len(feas_combined_partcoll))
        stiffrat_combined_nodalprop = np.zeros(len(feas_combined_nodalprop))
        stiffrat_combined_orient = np.zeros(len(feas_combined_orient))
        stiffrat_combined_inters = np.zeros(len(feas_combined_inters))
     
    # Compute Pareto Fronts
    pen_objs_combined = get_pen_objs(norm_objs_combined, feas_combined, conn_combined, stiffrat_combined)
    pen_objs_combined_partcoll = get_pen_objs(norm_objs_combined_partcoll, feas_combined_partcoll, conn_combined_partcoll, stiffrat_combined_partcoll)
    pen_objs_combined_nodalprop = get_pen_objs(norm_objs_combined_nodalprop, feas_combined_nodalprop, conn_combined_nodalprop, stiffrat_combined_nodalprop)
    pen_objs_combined_orient = get_pen_objs(norm_objs_combined_orient, feas_combined_orient, conn_combined_orient, stiffrat_combined_orient)
    pen_objs_combined_inters = get_pen_objs(norm_objs_combined_inters, feas_combined_inters, conn_combined_inters, stiffrat_combined_inters)
    
    pf_combined = compute_pareto_front(pen_objs_combined)
    pf_combined_partcoll = compute_pareto_front(pen_objs_combined_partcoll)
    pf_combined_nodalprop = compute_pareto_front(pen_objs_combined_nodalprop)
    pf_combined_orient = compute_pareto_front(pen_objs_combined_orient)
    pf_combined_inters = compute_pareto_front(pen_objs_combined_inters)
    
    data_comb = {}
    #data_comb['n_designs'] = len(norm_objs_combined)
    
    data_comb['pf'] = pf_combined
    data_comb['pf_partcoll'] = pf_combined_partcoll
    data_comb['pf_nodalprop'] = pf_combined_nodalprop
    data_comb['pf_orient'] = pf_combined_orient
    data_comb['pf_inters'] = pf_combined_inters
    
    return data_comb

def get_data_subset(data, dataset_inds, truss_prob):
    dataset = {}
    
    dataset['objs'] = data['objs'][dataset_inds,:]
    dataset['objs_partcoll'] = data['objs_partcoll'][dataset_inds,:]
    dataset['objs_nodalprop'] = data['objs_nodalprop'][dataset_inds,:]
    dataset['objs_orient'] = data['objs_orient'][dataset_inds,:]
    dataset['objs_inters'] = data['objs_inters'][dataset_inds,:]
    
    dataset['feas'] = data['feas'][dataset_inds]
    dataset['feas_partcoll'] = data['feas_partcoll'][dataset_inds]
    dataset['feas_nodalprop'] = data['feas_nodalprop'][dataset_inds]
    dataset['feas_orient'] = data['feas_orient'][dataset_inds]
    dataset['feas_inters'] = data['feas_inters'][dataset_inds]
    
    dataset['conn'] = data['conn'][dataset_inds]
    dataset['conn_partcoll'] = data['conn_partcoll'][dataset_inds]
    dataset['conn_nodalprop'] = data['conn_nodalprop'][dataset_inds]
    dataset['conn_orient'] = data['conn_orient'][dataset_inds]
    dataset['conn_inters'] = data['conn_inters'][dataset_inds]
    
    if truss_prob:
        dataset['stiffrat'] = data['stiffrat'][dataset_inds]
        dataset['stiffrat_partcoll'] = data['stiffrat_partcoll'][dataset_inds]
        dataset['stiffrat_nodalprop'] = data['stiffrat_nodalprop'][dataset_inds]
        dataset['stiffrat_orient'] = data['stiffrat_orient'][dataset_inds]
        dataset['stiffrat_inters'] = data['stiffrat_inters'][dataset_inds]
    
    return dataset
    

def get_obj_bounds_run(data_comb, truss_prob, rand_mode):
    #data_comb = get_combined_data(data_rand, data_moea, truss_prob, rand_mode)
    pf = data_comb['pf']
    pf_partcoll = data_comb['pf_partcoll']       
    pf_nodalprop = data_comb['pf_nodalprop']
    pf_orient = data_comb['pf_orient']
    pf_inters = data_comb['pf_inters']
    
    pf_objs1 = [x[0] for x in pf]
    pf_objs2 = [x[1] for x in pf]
    
    pf_partcoll_objs1 = [x[0] for x in pf_partcoll]
    pf_partcoll_objs2 = [x[1] for x in pf_partcoll]
    
    pf_nodalprop_objs1 = [x[0] for x in pf_nodalprop]
    pf_nodalprop_objs2 = [x[1] for x in pf_nodalprop]
    
    pf_orient_objs1 = [x[0] for x in pf_orient]
    pf_orient_objs2 = [x[1] for x in pf_orient]
    
    pf_inters_objs1 = [x[0] for x in pf_inters]
    pf_inters_objs2 = [x[1] for x in pf_inters]
    
    objs1_max = np.amax([np.amax(pf_objs1), np.amax(pf_partcoll_objs1), np.amax(pf_nodalprop_objs1), np.amax(pf_orient_objs1), np.amax(pf_inters_objs1)])
    objs1_min = np.amin([np.amin(pf_objs1), np.amin(pf_partcoll_objs1), np.amin(pf_nodalprop_objs1), np.amin(pf_orient_objs1), np.amin(pf_inters_objs1)])
    objs2_max = np.amax([np.amax(pf_objs2), np.amax(pf_partcoll_objs2), np.amax(pf_nodalprop_objs2), np.amax(pf_orient_objs2), np.amax(pf_inters_objs2)])
    objs2_min = np.amin([np.amin(pf_objs2), np.amin(pf_partcoll_objs2), np.amin(pf_nodalprop_objs2), np.amin(pf_orient_objs2), np.amin(pf_inters_objs2)])
    
    return [objs1_max, objs1_min, objs2_max, objs2_min]

def compute_I_hv(pf_h, pf):
    return compute_hv(pf_h) - compute_hv(pf)

def compute_indices(data_comb, truss_prob, obj_bounds_all): # to be called after combining data
    pf = data_comb['pf']
    pf_partcoll = data_comb['pf_partcoll']
    pf_nodalprop = data_comb['pf_nodalprop']
    pf_orient = data_comb['pf_orient']
    pf_inters = data_comb['pf_inters']

    #n_total = data_comb['n_designs']
    
    # Normalize pareto fronts
    pf_objs1_norm = [(x[0] - obj_bounds_all[1])/(obj_bounds_all[0] - obj_bounds_all[1]) for x in pf]
    pf_objs1_partcoll_norm = [(x[0] - obj_bounds_all[1])/(obj_bounds_all[0] - obj_bounds_all[1]) for x in pf_partcoll]
    pf_objs1_nodalprop_norm = [(x[0] - obj_bounds_all[1])/(obj_bounds_all[0] - obj_bounds_all[1]) for x in pf_nodalprop]
    pf_objs1_orient_norm = [(x[0] - obj_bounds_all[1])/(obj_bounds_all[0] - obj_bounds_all[1]) for x in pf_orient]
    pf_objs1_inters_norm = [(x[0] - obj_bounds_all[1])/(obj_bounds_all[0] - obj_bounds_all[1]) for x in pf_inters]
    
    pf_objs2_norm = [(x[1] - obj_bounds_all[3])/(obj_bounds_all[2] - obj_bounds_all[3]) for x in pf]
    pf_objs2_partcoll_norm = [(x[1] - obj_bounds_all[3])/(obj_bounds_all[2] - obj_bounds_all[3]) for x in pf_partcoll]
    pf_objs2_nodalprop_norm = [(x[1] - obj_bounds_all[3])/(obj_bounds_all[2] - obj_bounds_all[3]) for x in pf_nodalprop]
    pf_objs2_orient_norm = [(x[1] - obj_bounds_all[3])/(obj_bounds_all[2] - obj_bounds_all[3]) for x in pf_orient]
    pf_objs2_inters_norm = [(x[1] - obj_bounds_all[3])/(obj_bounds_all[2] - obj_bounds_all[3]) for x in pf_inters]
    
    pf_norm = np.column_stack((pf_objs1_norm, pf_objs2_norm))
    pf_partcoll_norm = np.column_stack((pf_objs1_partcoll_norm, pf_objs2_partcoll_norm))
    pf_nodalprop_norm = np.column_stack((pf_objs1_nodalprop_norm, pf_objs2_nodalprop_norm))
    pf_orient_norm = np.column_stack((pf_objs1_orient_norm, pf_objs2_orient_norm))
    pf_inters_norm = np.column_stack((pf_objs1_inters_norm, pf_objs2_inters_norm))
    
    # Compute indices based on objective minimization
    I_partcoll = compute_I_hv(pf_partcoll_norm, pf_norm) #compute_hv(pf_partcoll_norm) - compute_hv(pf_norm)
    I_nodalprop = compute_I_hv(pf_nodalprop_norm, pf_norm) #compute_hv(pf_nodalprop_norm) - compute_hv(pf_norm)
    I_orient = compute_I_hv(pf_orient_norm, pf_norm) #compute_hv(pf_orient_norm) - compute_hv(pf_norm)
    I_inters = compute_I_hv(pf_inters_norm, pf_norm) #compute_hv(pf_inters_norm) - compute_hv(pf_norm)
    
    return I_partcoll, I_nodalprop, I_orient, I_inters

############ OPERATION #############
### Compute the heuristic indices for each run
run_num = 0

dataset_mode = 3 # 1 - bootstrapping, 2 - leave one out, 3 - equal division

if dataset_mode == 1:
    n_datasets = 30
    n_des = 200
elif dataset_mode == 2:
    n_datasets = 300
    n_des = 299
else:
    n_datasets = 10
    n_des = 30

truss_problem = True
model_prob = 2 # 1 - fiber model, 2 - truss model, 3 - beam model
random_mode = 1 # 1 - only random data, 2 - random + MOEA data (always 1)

### Get all designs first
file_loc_rand_run = get_fileloc(truss_problem, model_prob, 1, run_num)
data_rand_run = read_csv(file_loc_rand_run, truss_problem)

### Compute Indices for each dataset
I_datasets = np.zeros((n_datasets, 4))
for i in range(n_datasets):
    
    if dataset_mode == 1:
        dataset_des_inds = np.random.randint(0, len(data_rand_run['objs']), n_des)
    elif dataset_mode == 2:
        dataset_des_inds_all = np.linspace(0, n_datasets-1, n_datasets)
        dataset_des_inds = np.delete(dataset_des_inds_all, i)
        dataset_des_inds = np.asarray(dataset_des_inds, dtype='int')
    else:
        dataset_des_inds = np.linspace(i*n_des, (i+1)*n_des-1, n_des)
        dataset_des_inds = np.asarray(dataset_des_inds, dtype='int')
        
    data_comb_run = get_data_subset(data_rand_run, dataset_des_inds, truss_problem)
    pfs_dataset_current = get_pfs_dataset(data_comb_run, truss_problem)
    obj_bounds_dataset_current = get_obj_bounds_run(pfs_dataset_current, truss_problem, random_mode)
    I_datasets[i,:] = compute_indices(pfs_dataset_current, truss_problem, obj_bounds_dataset_current)
    
### Empirical probability of positive index
prob_heurs = np.zeros((4))
for i in range(4):
    prob_heurs[i] = len([x for x in I_datasets[:,i] if x > 0])/n_datasets