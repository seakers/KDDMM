# -*- coding: utf-8 -*-
"""
Operator Index vs NFE for materials problems

@author: rosha
"""
from pygmo import hypervolume
import csv
import numpy as np
import matplotlib.pyplot as plt
from IPython.core.debugger import set_trace

### Useful functions
def find_closest_index(val, search_list):
    val_diff = np.array(search_list) - val
    closest_index = np.argmin(np.abs(val_diff))
    return closest_index

def find_last_index(val,search_list):
    if val in search_list:
        idx = len(search_list) - list(search_list[::-1]).index(val) - 1
    else:
        idx = 0
    return idx

def get_fileloc(truss_prob, model, run_mode, run_num):
    filepath_init = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\Operator Index Data\\' # for office system
    #filepath_init = 'C:\\Users\\rosha\\Documents\\SEAK Lab Github\\KD3M3\\result\\Operator Index Data\\' # for home system

    if run_mode == 1: 
        filepath_mode = 'Random\\'
        filename_init = 'random_operator_index_data_'
        filename_mode = ''
    elif run_mode == 2:
        filename_init = 'EpsilonMOEA_emoea__'
        filepath_mode = 'Epsilon MOEA\\'
        filename_mode = '_fullpop'

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
    
    filename = filename_init + filename_problem + filename_model + str(run_num) + filename_mode + '.csv'
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

def compute_pareto_front_constr(population_objs, population_constr_aggr):
    # population_constr_aggr -> sum of absolute constraint violations
    pop_size = len(population_objs)
    obj_num = 2

    domination_counter = [0] * pop_size

    for i in range(pop_size):
        for j in range(i+1, pop_size):
            # First check for aggregate constraint dominance
            if population_constr_aggr[i] < population_constr_aggr[j]:
                domination_counter[j] += 1
            elif population_constr_aggr[i] > population_constr_aggr[j]:
                domination_counter[i] += 1
            else:
                # For equal constraint satisfaction, check each objective for dominance
                dominate = [0] * obj_num
                for k in range(obj_num):
                    if population_objs[i][k] > population_objs[j][k]:
                        dominate[k] = 1
                    elif population_objs[i][k] < population_objs[j][k]:
                        dominate[k] = -1
                if -1 not in dominate and 1 in dominate:
                    domination_counter[i] += 1
                elif -1 in dominate and 1 not in dominate:
                    domination_counter[j] += 1
                
    pareto_solutions = []
    for i in range(len(domination_counter)):
        if domination_counter[i] == 0:
            pareto_solutions.append(population_objs[i])
    return pareto_solutions

def get_aggr_constr(feas_array, conn_array, stiffrat_array):
    feas_viol_array = np.subtract(1, feas_array)
    conn_viol_array = np.subtract(1, conn_array)
    aggr_constr1 = np.add(feas_viol_array, conn_viol_array)
    aggr_constr = np.add(aggr_constr1, stiffrat_array)
    return aggr_constr

def compute_hv(population):
    array_archs = np.zeros((len(population), 2))
    for i in range(len(population)):
        array_archs[i] = population[i]
    hv_object = hypervolume(array_archs)
    hv = hv_object.compute([1.1,1.1])/1.1**2
    return hv

def read_csv(fileloc_csv, truss_prob, run_mode, nfe_thresh):
    
    with open(fileloc_csv, newline='') as csvfile:
        data = [row for row in csv.reader(csvfile)]
        
        if run_mode == 2:
            nfes_dat = np.zeros(len(data)-1)
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
        
        buffer = 0
        if run_mode == 2:
            buffer = 1
        
        valid_count = 0
        for x in range(len(data)-1):
            if run_mode == 2:
                nfe_float = float(data[x+1][0])
            
            if truss_prob:
                data_float = list(map(float,data[x+1][buffer+1:buffer+6]))
                data_partcoll_float = list(map(float,data[x+1][buffer+7:buffer+12]))
                data_nodalprop_float = list(map(float,data[x+1][buffer+13:buffer+18]))
                data_orient_float = list(map(float,data[x+1][buffer+19:buffer+24]))
                data_inters_float = list(map(float,data[x+1][buffer+25:buffer+30]))
                
            else:
                data_float = list(map(float,data[x+1][buffer+1:buffer+5]))
                data_partcoll_float = list(map(float,data[x+1][buffer+6:buffer+10]))
                data_nodalprop_float = list(map(float,data[x+1][buffer+11:buffer+15]))
                data_orient_float = list(map(float,data[x+1][buffer+16:buffer+20]))
                data_inters_float = list(map(float,data[x+1][buffer+21:buffer+25]))
                
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
            
            if run_mode == 2:
                nfes_dat[valid_count] = nfe_float
                
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
                    stiffrat_dat[valid_count] = 10
                else:
                    stiffrat_dat[valid_count] = data_float[4]
                if (data_partcoll_float[4] > 50):
                    stiffrat_partcoll_dat[valid_count] = 10
                else:
                    stiffrat_partcoll_dat[valid_count] = data_partcoll_float[4]
                if (data_nodalprop_float[4] > 50):
                    stiffrat_nodalprop_dat[valid_count] = 10
                else:
                    stiffrat_nodalprop_dat[valid_count] = data_nodalprop_float[4]
                if (data_orient_float[4] > 50):
                    stiffrat_orient_dat[valid_count] = 10
                else:
                    stiffrat_orient_dat[valid_count] = data_orient_float[4]
                if (data_inters_float[4] > 50):
                    stiffrat_inters_dat[valid_count] = 10
                else:
                    stiffrat_inters_dat[valid_count] = data_inters_float[4]
                
            valid_count += 1
            
    if run_mode == 2:
        nfes = nfes_dat[:valid_count]
        
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
    
    # Sort by NFEs
    if run_mode == 2:
        sort_indices = np.argsort(nfes)
        
        nfes_sorted = nfes[sort_indices]
        
        norm_objs1_sorted = norm_objs1[sort_indices]
        norm_objs1_partcoll_sorted = norm_objs1_partcoll[sort_indices]
        norm_objs1_nodalprop_sorted = norm_objs1_nodalprop[sort_indices]
        norm_objs1_orient_sorted = norm_objs1_orient[sort_indices]
        norm_objs1_inters_sorted = norm_objs1_inters[sort_indices]
        
        norm_objs2_sorted = norm_objs2[sort_indices]
        norm_objs2_partcoll_sorted = norm_objs2_partcoll[sort_indices]
        norm_objs2_nodalprop_sorted = norm_objs2_nodalprop[sort_indices]
        norm_objs2_orient_sorted = norm_objs2_orient[sort_indices]
        norm_objs2_inters_sorted = norm_objs2_inters[sort_indices]
        
        feas_sorted = feas[sort_indices]
        feas_partcoll_sorted = feas_partcoll[sort_indices]
        feas_nodalprop_sorted = feas_nodalprop[sort_indices]
        feas_orient_sorted = feas_orient[sort_indices]
        feas_inters_sorted = feas_inters[sort_indices]
        
        conn_sorted = conn[sort_indices]
        conn_partcoll_sorted = conn_partcoll[sort_indices]
        conn_nodalprop_sorted = conn_nodalprop[sort_indices]
        conn_orient_sorted = conn_orient[sort_indices]
        conn_inters_sorted = conn_inters[sort_indices]
        
        if truss_prob:
             stiffrat_sorted = stiffrat[sort_indices]
             stiffrat_partcoll_sorted = stiffrat_partcoll[sort_indices]
             stiffrat_nodalprop_sorted = stiffrat_nodalprop[sort_indices]
             stiffrat_orient_sorted = stiffrat_orient[sort_indices]
             stiffrat_inters_sorted = stiffrat_inters[sort_indices]
    
        nfe_thresh_idx = find_last_index(nfe_thresh, nfes_sorted)
        
        norm_objs1 = norm_objs1_sorted[:nfe_thresh_idx]
        norm_objs1_partcoll = norm_objs1_partcoll_sorted[:nfe_thresh_idx]
        norm_objs1_nodalprop = norm_objs1_nodalprop_sorted[:nfe_thresh_idx]
        norm_objs1_orient = norm_objs1_orient_sorted[:nfe_thresh_idx]
        norm_objs1_inters = norm_objs1_inters_sorted[:nfe_thresh_idx]
        
        norm_objs2 = norm_objs2_sorted[:nfe_thresh_idx]
        norm_objs2_partcoll = norm_objs2_partcoll_sorted[:nfe_thresh_idx]
        norm_objs2_nodalprop = norm_objs2_nodalprop_sorted[:nfe_thresh_idx]
        norm_objs2_orient = norm_objs2_orient_sorted[:nfe_thresh_idx]
        norm_objs2_inters = norm_objs2_inters_sorted[:nfe_thresh_idx]
        
        feas = feas_sorted[:nfe_thresh_idx]
        feas_partcoll = feas_partcoll_sorted[:nfe_thresh_idx]
        feas_nodalprop = feas_nodalprop_sorted[:nfe_thresh_idx]
        feas_orient = feas_orient_sorted[:nfe_thresh_idx]
        feas_inters = feas_inters_sorted[:nfe_thresh_idx]
        
        conn = conn_sorted[:nfe_thresh_idx]
        conn_partcoll = conn_partcoll_sorted[:nfe_thresh_idx]
        conn_nodalprop = conn_nodalprop_sorted[:nfe_thresh_idx]
        conn_orient = conn_orient_sorted[:nfe_thresh_idx]
        conn_inters = conn_inters_sorted[:nfe_thresh_idx]
        
        if truss_prob:
            stiffrat = stiffrat_sorted[:nfe_thresh_idx]
            stiffrat_partcoll = stiffrat_partcoll_sorted[:nfe_thresh_idx]
            stiffrat_nodalprop = stiffrat_nodalprop_sorted[:nfe_thresh_idx]
            stiffrat_orient = stiffrat_orient_sorted[:nfe_thresh_idx]
            stiffrat_inters = stiffrat_inters_sorted[:nfe_thresh_idx]
        
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

def get_combined_data(data_rand, data_moea, truss_prob, rand_mode, hv_constr_only):
    
    if rand_mode == 1:
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
        norm_objs_combined = data_rand['objs'].tolist() + data_moea['objs'].tolist()
        norm_objs_combined_partcoll = data_rand['objs_partcoll'].tolist() + data_moea['objs_partcoll'].tolist()
        norm_objs_combined_nodalprop = data_rand['objs_nodalprop'].tolist() + data_moea['objs_nodalprop'].tolist()
        norm_objs_combined_orient = data_rand['objs_orient'].tolist() + data_moea['objs_orient'].tolist()
        norm_objs_combined_inters = data_rand['objs_inters'].tolist() + data_moea['objs_inters'].tolist()
        
        feas_combined = data_rand['feas'].tolist() + data_moea['feas'].tolist()
        feas_combined_partcoll = data_rand['feas_partcoll'].tolist() + data_moea['feas_partcoll'].tolist()
        feas_combined_nodalprop = data_rand['feas_nodalprop'].tolist() + data_moea['feas_nodalprop'].tolist()
        feas_combined_orient = data_rand['feas_orient'].tolist() + data_moea['feas_orient'].tolist()
        feas_combined_inters = data_rand['feas_inters'].tolist() + data_moea['feas_inters'].tolist()
        
        conn_combined = data_rand['conn'].tolist() + data_moea['conn'].tolist()
        conn_combined_partcoll = data_rand['conn_partcoll'].tolist() + data_moea['conn_partcoll'].tolist()
        conn_combined_nodalprop = data_rand['conn_nodalprop'].tolist() + data_moea['conn_nodalprop'].tolist()
        conn_combined_orient = data_rand['conn_orient'].tolist() + data_moea['conn_orient'].tolist()
        conn_combined_inters = data_rand['conn_inters'].tolist() + data_moea['conn_inters'].tolist()
        
        if truss_prob:
            stiffrat_combined = data_rand['stiffrat'].tolist() + data_moea['stiffrat'].tolist()
            stiffrat_combined_partcoll = data_rand['stiffrat_partcoll'].tolist() + data_moea['stiffrat_partcoll'].tolist()
            stiffrat_combined_nodalprop = data_rand['stiffrat_nodalprop'].tolist() + data_moea['stiffrat_nodalprop'].tolist()
            stiffrat_combined_orient = data_rand['stiffrat_orient'].tolist() + data_moea['stiffrat_orient'].tolist()
            stiffrat_combined_inters = data_rand['stiffrat_inters'].tolist() + data_moea['stiffrat_inters'].tolist()
        else:
            stiffrat_combined = np.zeros(len(feas_combined))
            stiffrat_combined_partcoll = np.zeros(len(feas_combined_partcoll))
            stiffrat_combined_nodalprop = np.zeros(len(feas_combined_nodalprop))
            stiffrat_combined_orient = np.zeros(len(feas_combined_orient))
            stiffrat_combined_inters = np.zeros(len(feas_combined_inters))
     
    # Compute Pareto Fronts
    if hv_constr_only:
        aggr_constr = get_aggr_constr(feas_combined, conn_combined, stiffrat_combined)
        aggr_constr_partcoll = get_aggr_constr(feas_combined_partcoll, conn_combined_partcoll, stiffrat_combined_partcoll)
        aggr_constr_nodalprop = get_aggr_constr(feas_combined_nodalprop, conn_combined_nodalprop, stiffrat_combined_nodalprop)
        aggr_constr_orient = get_aggr_constr(feas_combined_orient, conn_combined_orient, stiffrat_combined_orient)
        aggr_constr_inters = get_aggr_constr(feas_combined_inters, conn_combined_inters, stiffrat_combined_inters)
    
        pf_combined = compute_pareto_front_constr(norm_objs_combined, aggr_constr)
        pf_combined_partcoll = compute_pareto_front_constr(norm_objs_combined_partcoll, aggr_constr_partcoll)
        pf_combined_nodalprop = compute_pareto_front_constr(norm_objs_combined_nodalprop, aggr_constr_nodalprop)
        pf_combined_orient = compute_pareto_front_constr(norm_objs_combined_orient, aggr_constr_orient)
        pf_combined_inters = compute_pareto_front_constr(norm_objs_combined_inters, aggr_constr_inters)
        
    else:
        pf_combined = compute_pareto_front(norm_objs_combined)
        pf_combined_partcoll = compute_pareto_front(norm_objs_combined_partcoll)
        pf_combined_nodalprop = compute_pareto_front(norm_objs_combined_nodalprop)
        pf_combined_orient = compute_pareto_front(norm_objs_combined_orient)
        pf_combined_inters = compute_pareto_front(norm_objs_combined_inters)
    
    data_comb = {}
    data_comb['n_designs'] = len(norm_objs_combined)
    
    data_comb['pf'] = pf_combined
    data_comb['pf_partcoll'] = pf_combined_partcoll
    data_comb['pf_nodalprop'] = pf_combined_nodalprop
    data_comb['pf_orient'] = pf_combined_orient
    data_comb['pf_inters'] = pf_combined_inters
    
    if not hv_constr_only:
        data_comb['feas'] = feas_combined
        data_comb['feas_partcoll'] = feas_combined_partcoll
        data_comb['feas_nodalprop'] = feas_combined_nodalprop
        data_comb['feas_orient'] = feas_combined_orient
        data_comb['feas_inters'] = feas_combined_inters
    
        data_comb['conn'] = conn_combined
        data_comb['conn_partcoll'] = conn_combined_partcoll
        data_comb['conn_nodalprop'] = conn_combined_nodalprop
        data_comb['conn_orient'] = conn_combined_orient
        data_comb['conn_inters'] = conn_combined_inters
    
        if truss_prob:
            data_comb['stiffrat'] = stiffrat_combined
            data_comb['stiffrat_partcoll'] = stiffrat_combined_partcoll
            data_comb['stiffrat_nodalprop'] = stiffrat_combined_nodalprop
            data_comb['stiffrat_orient'] = stiffrat_combined_orient
            data_comb['stiffrat_inters'] = stiffrat_combined_inters
    
    return data_comb

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

def compute_I_hv2(pf_h, pf):
    return (compute_hv(pf_h)/compute_hv(pf)) - 1

def compute_I_constr(constrv_h, constrv, s_h, s):
    return (np.mean(constrv_h) - np.mean(constrv))*(s_h/s)

def compute_I_constr_alt(constrv_h, constrv):
    return (np.mean(constrv) - np.mean(constrv_h))/(np.max(constrv.tolist() + constrv_h.tolist()))

def compute_I_constr2(constrv_h, constrv):
    return 1 - (np.mean(constrv_h)/np.mean(constrv))
    
def compute_indices(data_comb, truss_prob, obj_bounds_all, hv_constr_only): # to be called after combining data
    pf = data_comb['pf']
    pf_partcoll = data_comb['pf_partcoll']
    pf_nodalprop = data_comb['pf_nodalprop']
    pf_orient = data_comb['pf_orient']
    pf_inters = data_comb['pf_inters']
    
    if not hv_constr_only:
        feas = data_comb['feas']
        feas_partcoll = data_comb['feas_partcoll']
        feas_nodalprop = data_comb['feas_nodalprop']
        feas_orient = data_comb['feas_orient']
        feas_inters = data_comb['feas_inters']
    
        conn = data_comb['conn']
        conn_partcoll = data_comb['conn_partcoll']
        conn_nodalprop = data_comb['conn_nodalprop']
        conn_orient = data_comb['conn_orient']
        conn_inters = data_comb['conn_inters']
        
        if truss_prob:
            stiffrat = data_comb['stiffrat']
            stiffrat_partcoll = data_comb['stiffrat_partcoll']
            stiffrat_nodalprop = data_comb['stiffrat_nodalprop']
            stiffrat_orient = data_comb['stiffrat_orient']
            stiffrat_inters = data_comb['stiffrat_inters']

    n_total = data_comb['n_designs']
    
    # Compute support of problem parameters
    if not hv_constr_only:
        s_pf = len(pf)/n_total
        s_pf_partcoll = len(pf_partcoll)/n_total
        s_pf_nodalprop = len(pf_nodalprop)/n_total
        s_pf_orient = len(pf_orient)/n_total
        s_pf_inters = len(pf_inters)/n_total
        
        s_feas = len([elem for elem in feas if elem==1])/n_total + 1e-5
        s_feas_partcoll = len([elem for elem in feas_partcoll if elem==1])/n_total + 1e-5
        s_feas_nodalprop = len([elem for elem in feas_nodalprop if elem==1])/n_total + 1e-5
        s_feas_orient = len([elem for elem in feas_orient if elem==1])/n_total + 1e-5
        s_feas_inters = len([elem for elem in feas_inters if elem==1])/n_total + 1e-5
        
        s_conn = len([elem for elem in conn if elem==1])/n_total + 1e-5
        s_conn_partcoll = len([elem for elem in conn_partcoll if elem==1])/n_total + 1e-5
        s_conn_nodalprop = len([elem for elem in conn_nodalprop if elem==1])/n_total + 1e-5
        s_conn_orient = len([elem for elem in conn_orient if elem==1])/n_total + 1e-5
        s_conn_inters = len([elem for elem in conn_inters if elem==1])/n_total + 1e-5
    
        if truss_prob:
            s_stiffrat = len([elem for elem in stiffrat if elem==0])/n_total + 1e-5
            s_stiffrat_partcoll = len([elem for elem in stiffrat_partcoll if elem==0])/n_total + 1e-5
            s_stiffrat_nodalprop = len([elem for elem in stiffrat_nodalprop if elem==0])/n_total + 1e-5
            s_stiffrat_orient = len([elem for elem in stiffrat_orient if elem==0])/n_total + 1e-5
            s_stiffrat_inters = len([elem for elem in stiffrat_inters if elem==0])/n_total + 1e-5
        
        s = [s_pf, s_feas, s_conn]
        s_partcoll = [s_pf_partcoll, s_feas_partcoll, s_conn_partcoll]
        s_nodalprop = [s_pf_nodalprop, s_feas_nodalprop, s_conn_nodalprop]
        s_orient = [s_pf_orient, s_feas_orient, s_conn_orient]
        s_inters = [s_pf_inters, s_feas_inters, s_conn_inters]
    
        if truss_prob:
            s.append(s_stiffrat)
            s_partcoll.append(s_stiffrat_partcoll)
            s_nodalprop.append(s_stiffrat_nodalprop)
            s_orient.append(s_stiffrat_orient)
            s_inters.append(s_stiffrat_inters)
        
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
        
    # Compute heuristic indices
    I_partcoll = 0
    I_nodalprop = 0
    I_orient = 0
    I_inters = 0    
    
    # Compute indices based on objective minimization
    I_partcoll = compute_I_hv(pf_partcoll_norm, pf_norm) #compute_hv(pf_partcoll_norm) - compute_hv(pf_norm)
    I_nodalprop = compute_I_hv(pf_nodalprop_norm, pf_norm) #compute_hv(pf_nodalprop_norm) - compute_hv(pf_norm)
    I_orient = compute_I_hv(pf_orient_norm, pf_norm) #compute_hv(pf_orient_norm) - compute_hv(pf_norm)
    I_inters = compute_I_hv(pf_inters_norm, pf_norm) #compute_hv(pf_inters_norm) - compute_hv(pf_norm)
    
    #I_partcoll = compute_I_hv2(pf_partcoll_norm, pf_norm) 
    #I_nodalprop = compute_I_hv2(pf_nodalprop_norm, pf_norm) 
    #I_orient = compute_I_hv2(pf_orient_norm, pf_norm) 
    #I_inters = compute_I_hv2(pf_inters_norm, pf_norm) 
    
    if not hv_constr_only:
        # Add indices based on feasibility satisfaction
        #I_partcoll += compute_I_constr(np.subtract(1, feas_partcoll), np.subtract(1, feas), s_feas_partcoll, s_feas) #(np.mean(np.subtract(1, feas)) - np.mean(np.subtract(1, feas_partcoll)))*(s_feas_partcoll/s_feas)
        #I_nodalprop += compute_I_constr(np.subtract(1, feas_nodalprop), np.subtract(1, feas), s_feas_nodalprop, s_feas) #(np.mean(np.subtract(1, feas)) - np.mean(np.subtract(1, feas_nodalprop)))*(s_feas_nodalprop/s_feas)
        #I_orient += compute_I_constr(np.subtract(1, feas_orient), np.subtract(1, feas), s_feas_orient, s_feas) #(np.mean(np.subtract(1, feas)) - np.mean(np.subtract(1, feas_orient)))*(s_feas_orient/s_feas) 
        #I_inters += compute_I_constr(np.subtract(1, feas_inters), np.subtract(1, feas), s_feas_inters, s_feas) #(np.mean(np.subtract(1, feas)) - np.mean(np.subtract(1, feas_inters)))*(s_feas_inters/s_feas) 
        
        I_partcoll += compute_I_constr2(np.subtract(1, feas_partcoll), np.subtract(1, feas))
        I_nodalprop += compute_I_constr2(np.subtract(1, feas_nodalprop), np.subtract(1, feas))
        I_orient += compute_I_constr2(np.subtract(1, feas_orient), np.subtract(1, feas))
        I_inters += compute_I_constr2(np.subtract(1, feas_inters), np.subtract(1, feas))
        
        # Add indices based on connectivity satisfaction
        #I_partcoll += compute_I_constr(np.subtract(1, conn_partcoll), np.subtract(1, conn), s_conn_partcoll, s_conn) #(np.mean(np.subtract(1, conn)) - np.mean(np.subtract(1, conn_partcoll)))*(s_conn_partcoll/s_conn) 
        #I_nodalprop += compute_I_constr(np.subtract(1, conn_nodalprop), np.subtract(1, conn), s_conn_nodalprop, s_conn) #(np.mean(np.subtract(1, conn)) - np.mean(np.subtract(1, conn_nodalprop)))*(s_conn_nodalprop/s_conn) 
        #I_orient += compute_I_constr(np.subtract(1, conn_orient), np.subtract(1, conn), s_conn_orient, s_conn) #(np.mean(np.subtract(1, conn)) - np.mean(np.subtract(1, conn_orient)))*(s_conn_orient/s_conn) 
        #I_inters += compute_I_constr(np.subtract(1, conn_inters), np.subtract(1, conn), s_conn_inters, s_conn) #(np.mean(np.subtract(1, conn)) - np.mean(np.subtract(1, conn_inters)))*(s_conn_inters/s_conn) 
        
        I_partcoll += compute_I_constr2(np.subtract(1, conn_partcoll), np.subtract(1, conn))
        I_nodalprop += compute_I_constr2(np.subtract(1, conn_nodalprop), np.subtract(1, conn))
        I_orient += compute_I_constr2(np.subtract(1, conn_orient), np.subtract(1, conn))
        I_inters += compute_I_constr2(np.subtract(1, conn_inters), np.subtract(1, conn))
    
        # If applicable, add indices based on stiffness ratio satisfaction
        if truss_prob:
            #I_partcoll += compute_I_constr(stiffrat_partcoll, stiffrat, s_stiffrat_partcoll, s_stiffrat) #(np.mean(stiffrat) - np.mean(stiffrat_partcoll))*(s_stiffrat_partcoll/s_stiffrat) 
            #I_nodalprop += compute_I_constr(stiffrat_nodalprop, stiffrat, s_stiffrat_nodalprop, s_stiffrat) #(np.mean(stiffrat) - np.mean(stiffrat_nodalprop))*(s_stiffrat_nodalprop/s_stiffrat) 
            #I_orient += compute_I_constr(stiffrat_orient, stiffrat, s_stiffrat_orient, s_stiffrat) #(np.mean(stiffrat) - np.mean(stiffrat_orient))*(s_stiffrat_orient/s_stiffrat) 
            #I_inters += compute_I_constr(stiffrat_inters, stiffrat, s_stiffrat_inters, s_stiffrat) #(np.mean(stiffrat) - np.mean(stiffrat_inters))*(s_stiffrat_inters/s_stiffrat) 
        
            I_partcoll += compute_I_constr2(stiffrat_partcoll, stiffrat) 
            I_nodalprop += compute_I_constr2(stiffrat_nodalprop, stiffrat) 
            I_orient += compute_I_constr2(stiffrat_orient, stiffrat) 
            I_inters += compute_I_constr2(stiffrat_inters, stiffrat) 
        
        if truss_prob:
            I_partcoll /= 4
            I_nodalprop /= 4
            I_orient /= 4
            I_inters /= 4
        else:
            I_partcoll /= 3
            I_nodalprop /= 3
            I_orient /= 3
            I_inters /= 3
    
    return I_partcoll, I_nodalprop, I_orient, I_inters

def compute_indices_allruns(n_runs, truss_prob, model_prob, nfe_threshold, hv_constr_only):
    I_pc_allruns = np.zeros(n_runs)
    I_np_allruns = np.zeros(n_runs)
    I_orient_allruns = np.zeros(n_runs)
    I_inters_allruns = np.zeros(n_runs)
    
    ### Compute Pareto Fronts for all runs first
    objs1_max_allruns = np.zeros(n_runs)
    objs1_min_allruns = np.zeros(n_runs)
    objs2_max_allruns = np.zeros(n_runs)
    objs2_min_allruns = np.zeros(n_runs)
    data_comb_allruns = {}
    for i in range(n_runs):
        file_loc_rand_i = get_fileloc(truss_prob, model_prob, 1, i)
        data_rand_i = read_csv(file_loc_rand_i, truss_prob, 1, 0)
        
        file_loc_moea_i = get_fileloc(truss_prob, model_prob, 2, i)
        data_moea_i = read_csv(file_loc_moea_i, truss_prob, 2, nfe_threshold)

        data_comb_i = get_combined_data(data_rand_i, data_moea_i, truss_prob, 2, hv_constr_only)
        
        data_comb_allruns['run'+str(i)] = data_comb_i
        
        obj_bounds_i = get_obj_bounds_run(data_comb_i, truss_prob, 2)
        objs1_max_allruns[i] = obj_bounds_i[0]
        objs1_min_allruns[i] = obj_bounds_i[1]
        objs2_max_allruns[i] = obj_bounds_i[2]
        objs2_min_allruns[i] = obj_bounds_i[3]
    
    ### Compute overall objective bounds 
    objs1_max_overall = np.amax(objs1_max_allruns)
    objs1_min_overall = np.amin(objs1_min_allruns)
    objs2_max_overall = np.amax(objs2_max_allruns)
    objs2_min_overall = np.amin(objs2_min_allruns)
    obj_bounds_overall = [objs1_max_overall, objs1_min_overall, objs2_max_overall, objs2_min_overall] 

    ### Compute indices
    for i in range(n_runs):
        data_comb_i = data_comb_allruns['run'+str(i)]
        I_pci, I_npi, I_orienti, I_intersi = compute_indices(data_comb_i, truss_prob, obj_bounds_overall, hv_constr_only)
        I_pc_allruns[i] = I_pci
        I_np_allruns[i] = I_npi
        I_orient_allruns[i] = I_orienti
        I_inters_allruns[i] = I_intersi
    
    return I_pc_allruns, I_np_allruns, I_orient_allruns, I_inters_allruns

### Operation
num_runs = 10
nfe_max = 6000
n_datapoints = 10

truss_problem = False
model_problem = 2 # 1 - fiber model, 2 - truss model, 3 - beam model
only_hv_constr = True # Use only HV of constrained Pareto Fronts for indices computation

nfes_array_f = np.linspace(0, nfe_max, n_datapoints+1)
nfes_array = [int(x) for x in nfes_array_f]

I_partcoll_allnfes = []
I_nodalprop_allnfes = []
I_orient_allnfes = []
I_inters_allnfes = []

for nfe in nfes_array:
    I_pc_nfe, I_np_nfe, I_orient_nfe, I_inters_nfe = compute_indices_allruns(num_runs, truss_problem, model_problem, nfe, only_hv_constr)
    I_partcoll_allnfes.append(I_pc_nfe)
    I_nodalprop_allnfes.append(I_np_nfe)
    I_orient_allnfes.append(I_orient_nfe)
    I_inters_allnfes.append(I_inters_nfe)
    
### Plotting
#fig, axs = plt.subplots(2, 2)

nfe_idx_ticks = np.linspace(1, len(nfes_array), len(nfes_array))

#plt.setp(axs, xticks=nfe_idx_ticks, xticklabels=nfes_array, ylabel='Heuristic Index')

#plt.sca(axs[0,0])
plt.figure()
plt.boxplot(I_partcoll_allnfes)
plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
plt.xlabel('NFE')
plt.ylabel('Heuristic Index', fontsize=8)
plt.title('Partial Collapsibility', fontsize=8)

#plt.sca(axs[0,1])
plt.figure()
plt.boxplot(I_nodalprop_allnfes)
plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
plt.xlabel('NFE')
plt.ylabel('Heuristic Index', fontsize=8)
plt.title('Nodal Properties', fontsize=8)

#plt.sca(axs[1,0])
plt.figure()
plt.boxplot(I_orient_allnfes)
plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
plt.xlabel('NFE')
plt.ylabel('Heuristic Index', fontsize=8)
plt.title('Orientation', fontsize=8)

#plt.sca(axs[1,1])
plt.figure()
plt.boxplot(I_inters_allnfes)
plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
plt.xlabel('NFE')
plt.ylabel('Heuristic Index', fontsize=8)
plt.title('Intersection', fontsize=8)

#fig.tight_layout(pad=3.0)

#plt.subplot(2,2,1)
#plt.boxplot(I_partcoll_allnfes)
#plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
#plt.ylabel('Operator Index', fontsize=8)
#plt.title('Partial Collapsibility', fontsize=8)

#plt.subplot(2,2,2)
#plt.boxplot(I_nodalprop_allnfes)
#plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
#plt.ylabel('Operator Index', fontsize=8)
#plt.title('Nodal Properties', fontsize=8)

#plt.subplot(2,2,3)
#plt.boxplot(I_orient_allnfes)
#plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
#plt.ylabel('Operator Index', fontsize=8)
#plt.title('Orientation', fontsize=8)

#plt.subplot(2,2,4)
#plt.boxplot(I_inters_allnfes)
#plt.xticks(nfe_idx_ticks, nfes_array, fontsize=8)
#plt.ylabel('Operator Index', fontsize=8)
#plt.title('Intersection', fontsize=8)

#plt.show()