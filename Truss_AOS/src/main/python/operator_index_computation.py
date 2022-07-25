# -*- coding: utf-8 -*-
"""
Compute Operator Index for materials problems

@author: rosha
"""
from pygmo import hypervolume
import csv
import math
import numpy as np
from IPython.core.debugger import set_trace

### Useful functions
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

def compute_hv(population):
    array_archs = np.zeros((len(population), 2))
    for i in range(len(population)):
        array_archs[i] = population[i]
    hv_object = hypervolume(array_archs)
    hv = hv_object.compute([1.1,1.1])/1.1**2
    return hv

def read_csv_and_compute_index(truss_problem, model, run_num):
    
    #filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\' # for office system
    filepath_init = 'C:\\Users\\rosha\\Documents\\SEAK Lab Github\\KD3M3\\result\\Operator Index Data\\' # for home system
    
    filename_init = 'random_operator_index_data_'
    
    filepath_problem = 'Artery Problem\\'
    filename_problem = 'artery_'
    if truss_problem:
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
    filepath = filepath_init + filepath_problem + filepath_model + filename
    
    with open(filepath,newline='') as csvfile:
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
        
        if truss_problem:
            stiffrat_dat = np.zeros(len(data)-1)
            stiffrat_partcoll_dat = np.zeros(len(data)-1)
            stiffrat_nodalprop_dat = np.zeros(len(data)-1)
            stiffrat_orient_dat = np.zeros(len(data)-1)
            stiffrat_inters_dat = np.zeros(len(data)-1)
        
        valid_count = 0
        for x in range(len(data)-1):
            if truss_problem:
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
            
            if truss_problem:
                stiffrat_dat[valid_count] = data_float[4]
                stiffrat_partcoll_dat[valid_count] = data_partcoll_float[4]
                stiffrat_nodalprop_dat[valid_count] = data_nodalprop_float[4]
                stiffrat_orient_dat[valid_count] = data_orient_float[4]
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
        
    if truss_problem:
        stiffrat = stiffrat_dat[:valid_count]
        stiffrat_partcoll = stiffrat_partcoll_dat[:valid_count]
        stiffrat_nodalprop = stiffrat_nodalprop_dat[:valid_count]
        stiffrat_orient = stiffrat_orient_dat[:valid_count]
        stiffrat_inters = stiffrat_inters_dat[:valid_count]
        
    ## Compute heuristic indices
    I_partcoll = 0
    I_nodalprop = 0
    I_orient = 0
    I_inters = 0
         
    # Compute Pareto Fronts
    norm_objs = np.column_stack((norm_objs1, norm_objs2))
    norm_objs_partcoll = np.column_stack((norm_objs1_partcoll, norm_objs2_partcoll))
    norm_objs_nodalprop = np.column_stack((norm_objs1_nodalprop, norm_objs2_nodalprop))
    norm_objs_orient = np.column_stack((norm_objs1_orient, norm_objs2_orient))
    norm_objs_inters = np.column_stack((norm_objs1_inters, norm_objs2_inters))
    
    pf = compute_pareto_front(norm_objs)
    pf_partcoll = compute_pareto_front(norm_objs_partcoll)
    pf_nodalprop = compute_pareto_front(norm_objs_nodalprop)
    pf_orient = compute_pareto_front(norm_objs_orient)
    pf_inters = compute_pareto_front(norm_objs_inters)
    
    ## Compute support of problem parameters
    s_pf = len(pf)/(len(data)-1)
    s_pf_partcoll = len(pf_partcoll)/(len(data)-1)
    s_pf_nodalprop = len(pf_nodalprop)/(len(data)-1)
    s_pf_orient = len(pf_orient)/(len(data)-1)
    s_pf_inters = len(pf_inters)/(len(data)-1)
    
    s_feas = len([elem for elem in feas if elem==1])/(len(data)-1) + 1e-5
    s_feas_partcoll = len([elem for elem in feas_partcoll if elem==1])/(len(data)-1) + 1e-5
    s_feas_nodalprop = len([elem for elem in feas_nodalprop if elem==1])/(len(data)-1) + 1e-5
    s_feas_orient = len([elem for elem in feas_orient if elem==1])/(len(data)-1) + 1e-5
    s_feas_inters = len([elem for elem in feas_inters if elem==1])/(len(data)-1) + 1e-5
    
    s_conn = len([elem for elem in conn if elem==1])/(len(data)-1) + 1e-5
    s_conn_partcoll = len([elem for elem in conn_partcoll if elem==1])/(len(data)-1) + 1e-5
    s_conn_nodalprop = len([elem for elem in conn_nodalprop if elem==1])/(len(data)-1) + 1e-5
    s_conn_orient = len([elem for elem in conn_orient if elem==1])/(len(data)-1) + 1e-5
    s_conn_inters = len([elem for elem in conn_inters if elem==1])/(len(data)-1) + 1e-5
    
    if truss_problem:
        s_stiffrat = len([elem for elem in stiffrat if elem==0])/(len(data)-1) + 1e-5
        s_stiffrat_partcoll = len([elem for elem in stiffrat_partcoll if elem==0])/(len(data)-1) + 1e-5
        s_stiffrat_nodalprop = len([elem for elem in stiffrat_nodalprop if elem==0])/(len(data)-1) + 1e-5
        s_stiffrat_orient = len([elem for elem in stiffrat_orient if elem==0])/(len(data)-1) + 1e-5
        s_stiffrat_inters = len([elem for elem in stiffrat_inters if elem==0])/(len(data)-1) + 1e-5
    
    # Compute indices based on objective minimization
    I_partcoll = compute_hv(pf_partcoll)*(-math.log10(s_pf_partcoll)) - compute_hv(pf)*(-math.log10(s_pf))
    I_nodalprop = compute_hv(pf_nodalprop)*(-math.log10(s_pf_nodalprop)) - compute_hv(pf)*(-math.log10(s_pf))
    I_orient = compute_hv(pf_orient)*(-math.log10(s_pf_orient)) - compute_hv(pf)*(-math.log10(s_pf))
    I_inters = compute_hv(pf_inters)*(-math.log10(s_pf_inters)) - compute_hv(pf)*(-math.log10(s_pf))
    
    # Add indices based on feasibility satisfaction
    I_partcoll += np.mean(np.subtract(1, feas))*(-math.log10(s_feas)) - np.mean(np.subtract(1, feas_partcoll))*(-math.log10(s_feas_partcoll))
    I_nodalprop += np.mean(np.subtract(1, feas))*(-math.log10(s_feas)) - np.mean(np.subtract(1, feas_nodalprop))*(-math.log10(s_feas_nodalprop))
    I_orient += np.mean(np.subtract(1, feas))*(-math.log10(s_feas)) - np.mean(np.subtract(1, feas_orient))*(-math.log10(s_feas_orient)) 
    I_inters += np.mean(np.subtract(1, feas))*(-math.log10(s_feas)) - np.mean(np.subtract(1, feas_inters))*(-math.log10(s_feas_inters)) 
    
    # Add indices based on connectivity satisfaction
    I_partcoll += np.mean(np.subtract(1, conn))*(-math.log10(s_conn)) - np.mean(np.subtract(1, conn_partcoll))*(-math.log10(s_conn_partcoll)) 
    I_nodalprop += np.mean(np.subtract(1, conn))*(-math.log10(s_conn)) - np.mean(np.subtract(1, conn_nodalprop))*(-math.log10(s_conn_nodalprop)) 
    I_orient += np.mean(np.subtract(1, conn))*(-math.log10(s_conn)) - np.mean(np.subtract(1, conn_orient))*(-math.log10(s_conn_orient)) 
    I_inters += np.mean(np.subtract(1, conn))*(-math.log10(s_conn)) - np.mean(np.subtract(1, conn_inters))*(-math.log10(s_conn_inters)) 
    
    # If applicable, add indices based on stiffness ratio satisfaction
    if truss_problem:
        I_partcoll += np.mean(stiffrat)*(-math.log10(s_stiffrat)) - np.mean(stiffrat_partcoll)*(-math.log10(s_stiffrat_partcoll)) 
        I_nodalprop += np.mean(stiffrat)*(-math.log10(s_stiffrat)) - np.mean(stiffrat_nodalprop)*(-math.log10(s_stiffrat_nodalprop)) 
        I_orient += np.mean(stiffrat)*(-math.log10(s_stiffrat)) - np.mean(stiffrat_orient)*(-math.log10(s_stiffrat_orient)) 
        I_inters += np.mean(stiffrat)*(-math.log10(s_stiffrat)) - np.mean(stiffrat_inters)*(-math.log10(s_stiffrat_inters)) 
        
    if truss_problem:
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

### Compute the heuristic indices for each run
n_runs = 10

I_pc = 0
I_np = 0
I_orient = 0
I_inters = 0

truss_prob = True
model_prob = 2 # 1 - fiber model, 2 - truss model, 3 - beam model

for i in range(n_runs):
    I_pci, I_npi, I_orienti, I_intersi = read_csv_and_compute_index(truss_prob, model_prob, i)
    I_pc += I_pci
    I_np += I_npi
    I_orient += I_orienti
    I_inters += I_intersi
    
I_pc /= n_runs
I_np /= n_runs
I_orient /= n_runs
I_inters /= n_runs

print('I_partcoll = ' + str(I_pc))
print('I_nodalprop = ' + str(I_np))
print('I_orient = ' + str(I_orient))
print('I_inters = ' + str(I_inters))
