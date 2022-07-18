# -*- coding: utf-8 -*-
"""
Hypervolume computation for both materials problems

@author: roshan94
"""
from pygmo import hypervolume
import csv
import statistics
import math
import numpy as np
import operator as op
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
from functools import reduce
from itertools import combinations
import matplotlib.pyplot as plt
import timeit

#### Useful functions and parameter defintions 
def nchoosek(n,k):
    k = min(k, n-k)
    num = reduce(op.mul, range(n,n-k,-1), 1)
    den = reduce(op.mul, range(1,k+1), 1)
    return num/den    
    
def get_true_objectives(true_obj1_array, true_obj2_array, index):
    return true_obj1_array[index], true_obj2_array[index]

def find_last_index(val,search_list):
    if val in search_list:
        idx = len(search_list) - search_list[::-1].index(val) - 1
    else:
        idx = 0
    return idx

def find_closest_index(val,search_list):
    val_diff = np.array(search_list) - val
    closest_index = np.argmin(np.abs(val_diff))
    return closest_index

def compute_pareto_front(population_objs, population_constr_aggr):
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

def compute_hv(population):
    array_archs = np.zeros((len(population), 2))
    for i in range(len(population)):
        array_archs[i] = population[i]
    hv_object = hypervolume(array_archs)
    hv = hv_object.compute([1.1,1.1])/1.1**2
    return hv

def get_array_element(array, index):
    return array[index]

#### Define NFE array for hypervolume computation (based on number of evaluations in optimization runs)
start_time = timeit.default_timer()

np.set_printoptions(threshold=np.inf)

# Create array of NFE values at which to compute hypervolume (assumes max function evaluations is 6000)
n_iter_total = 100 # Total number of points in NFE array (1 more than input value to incorporate 0)
n_iter_init = 80 # Number of initial points in NFE array separated by 50 (the rest after that are separated by 100)
nfe_array = np.zeros(n_iter_total+1)
for i in range(n_iter_init):
    nfe_array[i] = 50*i
    
for i in range(n_iter_total - n_iter_init + 1):
    nfe_array[n_iter_init+i] = 50*n_iter_init + 100*i

#### Determine csv filepath from given case type for one of the materials problems
def get_csv_filepath_material(partcoll_constrained, nodprop_constrained, orient_constrained, inters_constrained, artery_problem, model_choice, cred_assign, run_number):
    # partcoll_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # nodprop_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # orient_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # inters_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # artery_problem = True if artery problem data is to be read, False if truss problem data is to be read
    # model_choice = 1 - fibre stiffness, 2 - truss stiffness, 3 - ANSYS APDL beam model
    
    filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\' # for office system
    #filepath = 'C:\\Users\\rosha\\Documents\\SEAK Lab Github\\KD3M3\\result\\' # for home system
    methods = ['Int Pen','AOS','Bias Init','ACH']
    heurs_list = ['PartColl','NodalProp','Orient','Inters']
    heur_abbrvs_list = ['p','n','o','i']
    heur_bools = np.vstack((partcoll_constrained, nodprop_constrained, orient_constrained, inters_constrained))
    aos_bools = [x[1] for x in heur_bools]
    
    if artery_problem:
        filepath_prob = 'Artery Problem\\'
        filename_prob = '_artery'
    else:
        filepath_prob = 'Truss Problem\\'
        filename_prob = '_prob2'
        
    filepath_constrad = 'Constant Radii\\'
        
    filepath_cred = ''
    if (any(aos_bools)):
        filename = 'AOSMOEA_emoea_'
        filepath_cred = 'offspring parent dominance\\'
        if cred_assign == 1:
            filepath_cred = 'set improvement dominance\\'
        elif cred_assign == 2:
            filepath_cred = 'set contribution dominance\\'
    else:
        filename = 'EpsilonMOEA_emoea_'
        
    filepath2 = ''
    filename2 = ''
    constr_count = 0
    for i in range(len(heur_bools[0])):
        constraints = methods[i] + ' - '
        constraints_abbrv = ''
        heur_count = 0
        for j in range(len(heur_bools)):
            if heur_bools[j][i]:
                constraints = constraints + heurs_list[j]
                constraints_abbrv = constraints_abbrv + heur_abbrvs_list[j]
            else:
                heur_count += 1
            
        if heur_count < len(heur_bools):
            filepath2 = filepath2 + constraints + '\\'
            filename2 = filename2 + constraints_abbrv + 'con' + str(i) + '_'
        else:
            constr_count += 1
            
    filepath_moea = ''
    if (constr_count == len(heur_bools[0])):
        filepath_moea = 'Epsilon MOEA\\'
        
    if model_choice == 1:
        filepath_model = 'Fibre Model\\'
        filename_model = '_fibre_fullpop.csv'
    elif model_choice == 2:
        filepath_model = 'Truss Model\\'
        filename_model = '_truss_fullpop.csv'
    elif model_choice == 3:
        filepath_model = 'Beam Model\\'
        filename_model = '_beam_fullpop.csv'
        
    return filepath + filepath_prob + filepath_constrad + filepath_model + filepath2 + filepath_moea + filepath_cred + filename + str(run_number) + filename2 + filename_prob + filename_model

#### Extract Pareto Front and normalization constants data from csv file
def extract_data_from_csv(csv_filepath, artery_problem, intpen_constr_heur, sidenum):
    # intpen_constr_heur = [intpen_constr_partcoll, intpen_constr_nodalprop, intpen_constr_orient, intpen_constr_inters] boolean array
    n_total_members = int(nchoosek(sidenum**2,2))
    with open(csv_filepath,newline='') as csvfile:
        data = [row for row in csv.reader(csvfile)]
                
        num_func_evals_dat = np.zeros(len(data)-1)
        pen_obj1_dat = np.zeros(len(data)-1)
        pen_obj2_dat = np.zeros(len(data)-1)
        true_obj1_dat = np.zeros(len(data)-1)
        true_obj2_dat = np.zeros(len(data)-1)
        
        feas_scores_dat = np.zeros(len(data)-1)
        conn_scores_dat = np.zeros(len(data)-1)
        if (not artery_problem):
            stiffrat_vals_dat = np.zeros(len(data)-1)
        partcoll_scores_dat = np.zeros(len(data)-1)
        nodprop_scores_dat = np.zeros(len(data)-1)
        orient_scores_dat = np.zeros(len(data)-1)
        inters_scores_dat = np.zeros(len(data)-1)
        
        valid_count = 0
        for x in range(len(data)-1):
            data_float = list(map(float,data[x+1][1:]))
            if (any(np.isnan(np.array(data_float))) or any(np.isinf(np.array(data_float)))):
                continue
            
            num_func_evals_dat[valid_count] = int(data[x+1][1])
            pen_obj1_dat[valid_count] = float(data[x+1][2])
            pen_obj2_dat[valid_count] = float(data[x+1][3])
            true_obj1_dat[valid_count] = float(data[x+1][4])
            true_obj2_dat[valid_count] = float(data[x+1][5])
            
            feas_scores_dat[valid_count] = float(data[x+1][6])
            conn_scores_dat[valid_count] = float(data[x+1][7])
            if (not artery_problem):
                stiffrat_vals_dat[valid_count] = float(data[x+1][8])
            
            partcoll_scores_dat[valid_count] = float(data[x+1][-4])
            nodprop_scores_dat[valid_count] = float(data[x+1][-3])
            orient_scores_dat[valid_count] = float(data[x+1][-2])
            inters_scores_dat[valid_count] = float(data[x+1][-1])
    
            valid_count += 1
            
    #designs = designs_dat[:valid_count]
    num_func_evals = num_func_evals_dat[:valid_count]
    pen_obj1 = pen_obj1_dat[:valid_count]
    pen_obj2 = pen_obj2_dat[:valid_count]
    true_obj1 = true_obj1_dat[:valid_count]
    true_obj2 = true_obj2_dat[:valid_count]
    
    feas_scores =  feas_scores_dat[:valid_count]
    conn_scores = conn_scores_dat[:valid_count]
    if (not artery_problem):
        stiffrat_vals = stiffrat_vals_dat[:valid_count]
    partcoll_scores = partcoll_scores_dat[:valid_count]
    nodprop_scores = nodprop_scores_dat[:valid_count]
    orient_scores = orient_scores_dat[:valid_count]
    inters_scores = inters_scores_dat[:valid_count]
            
    #print('pen_obj1')
    #print(pen_obj1)
    #print('\n')
    #print('pen_obj2')
    #print(pen_obj2)
    #print('\n')
    #print('true_obj1')
    #print(true_obj1)
    #print('\n')
    #print('true_obj2')
    #print(true_obj2)
    #print('\n')
    
    ## Sort num_fun_evals (and obj1 & obj2, feas and stab scores) in ascending order
    n_func_evals = num_func_evals
    sort_indices = np.argsort(n_func_evals)
    pen_obj1_sorted = list(pen_obj1[sort_indices])
    pen_obj2_sorted = list(pen_obj2[sort_indices])
    true_obj1_sorted = list(true_obj1[sort_indices])
    true_obj2_sorted = list(true_obj2[sort_indices])
    
    feas_scores_sorted = list(feas_scores[sort_indices])
    conn_scores_sorted = list(conn_scores[sort_indices])
    if (not artery_problem):
        stiffrat_vals_sorted = list(stiffrat_vals[sort_indices])
    partcoll_scores_sorted = list(partcoll_scores[sort_indices])
    nodprop_scores_sorted = list(nodprop_scores[sort_indices])
    orient_scores_sorted = list(orient_scores[sort_indices])
    inters_scores_sorted = list(inters_scores[sort_indices])
    
    #designs_sorted = []
    #for i in range(len(sort_indices)):
        #designs_sorted.append(designs[sort_indices[i]])
    
    heur_scores_sorted = np.vstack((partcoll_scores_sorted, nodprop_scores_sorted, orient_scores_sorted, inters_scores_sorted))
    
    ## Compute only constraint optimized (normalized) objectives (used for the interior penalty cases)
    heur_weight = 1
    heur_pen = np.zeros(len(feas_scores_sorted))
    if (any(intpen_constr_heur)):
        heur_index_array = np.arange(len(intpen_constr_heur))
        heur_index_used = [v for i, v in enumerate(heur_index_array) if intpen_constr_heur[i] == True]
     
        for idx in heur_index_used:
            heur_scores_idx = heur_scores_sorted[idx]
            heur_pen_idx = [math.log10(abs(x))/16 for x in heur_scores_idx]
            heur_pen = np.add(heur_pen,np.divide(heur_pen_idx,len(heur_index_used)))
            
    weighted_heur_pen = [k*heur_weight for k in heur_pen]
    
    pen_obj1_constr_sorted = list(np.add(pen_obj1_sorted,weighted_heur_pen))
    pen_obj2_constr_sorted = list(np.add(pen_obj2_sorted,weighted_heur_pen))
    
    ## Determine normalizing objective scores and compute pareto fronts for optimized (normalized) and true objectives as well as for true objectives of only feasible designs 
    nfe_list_sorted = list(n_func_evals[sort_indices])
    
    #max_func_evals = nfe_list_sorted[-1]
    max_func_evals = 6000 # some runs for some cases run upto 5001 evaluations, which causes hv array length issues

    fullsat = []   
    pareto_front_dict = {}
    #pareto_front_feas_dict = {}
    #pareto_front_conn_dict = {}
    #pareto_front_stiffrat_dict = {}
    #pareto_front_partcoll_dict = {}
    #pareto_front_nodprop_dict = {}
    #pareto_front_orient_dict = {}
    #pareto_front_designs_dict = {}
    pareto_front_true_dict = {}
    pareto_front_truesat_dict = {}
    pf_normalize_max_obj1 = []
    pf_normalize_min_obj1 = []
    pf_normalize_max_obj2 = []
    pf_normalize_min_obj2 = []
    pf_true_normalize_max_obj1 = []
    pf_true_normalize_min_obj1 = []
    pf_true_normalize_max_obj2 = []
    pf_true_normalize_min_obj2 = []
    pf_truesat_normalize_max_obj1 = []
    pf_truesat_normalize_min_obj1 = []
    pf_truesat_normalize_max_obj2 = []
    pf_truesat_normalize_min_obj2 = []
    
    num_fully_satisfying = np.zeros(len(nfe_array))
    
    count = 0
    pop_size = int(find_last_index(0, nfe_list_sorted))
    nfe_jump_recorded = False
    jump_nfe = 0

    for i in range(len(nfe_array)):
        #print('iter = ' + str(i))
        nfe_val = nfe_array[i]
    
        if (nfe_list_sorted[0] == 0):
            if (nfe_val <= 100): # population size set at 100 in java code, but maybe different here due to NaNs
                nfe_index_current = pop_size+1
                nfe_index_previous = 0 #NEW
            else:
                nfe_index_current = find_closest_index(nfe_val, nfe_list_sorted)
                nfe_index_previous = find_closest_index(nfe_array[i-1], nfe_list_sorted) #NEW         
        else:
            if (nfe_val <= nfe_list_sorted[0]):
                nfe_index_previous = 0 #NEW
                nfe_index_current = find_closest_index(0, nfe_list_sorted)
            else:
                nfe_index_current = find_closest_index(nfe_val, nfe_list_sorted)
                nfe_index_previous = find_closest_index(nfe_array[i-1], nfe_list_sorted) #NEW
        
        nfe_array_current = nfe_list_sorted[nfe_index_previous:nfe_index_current] #NEW
        #nfe_array_current = nfe_list_sorted[:nfe_index_current]
        current_population = []
        for j in range(len(nfe_array_current)):
            current_population.append([pen_obj1_constr_sorted[nfe_index_previous+j], pen_obj2_constr_sorted[nfe_index_previous+j]])
            
        current_population_unique = np.unique(current_population, axis=0).tolist()
        
        current_population_constr_aggr = []
        current_population_truesat = []
        for design in current_population_unique:
            des_index = pen_obj1_constr_sorted.index(design[0])
            des_constr = []
            des_constr.append(1 - get_array_element(feas_scores_sorted, des_index))
            des_constr.append(1 - get_array_element(conn_scores_sorted, des_index))
            if (not artery_problem):
                des_constr.append(get_array_element(stiffrat_vals_sorted, des_index))
            
            constr_aggr = np.sum(np.abs(des_constr))
            current_population_constr_aggr.append(constr_aggr)
            if (constr_aggr == 0):
                true_obj1, true_obj2 = get_true_objectives(true_obj1_sorted, true_obj2_sorted, des_index)
                current_population_truesat.append([-true_obj1, true_obj2])
        
        # if (not current_population_truesat):
        #     # if (i == 0):
        #     #     num_fully_satisfying[i] = 0
        #     # else:
        #     #     num_fully_satisfying[i] = num_fully_satisfying[i-1] + 0
        #     #num_fully_satisfying[i] = 0
        # else:
        #     current_population_truesat_unique = np.unique(current_population_truesat, axis = 0)
        #     #for k2 in range(len(current_population_truesat_unique)):
        #         #fullsat.append(current_population_truesat_unique[k2])
            
        #     #fullsat = np.unique(fullsat, axis = 0).tolist()
            
        #     # if (i == 0):
        #     #     num_fully_satisfying[i] = len(current_population_truesat_unique)
        #     # else:
        #     #     num_fully_satisfying[i] = num_fully_satisfying[i-1] + len(current_population_truesat_unique)
        #     ##num_fully_satisfying[i] = len(current_population_truesat_unique)
        #     #num_fully_satisfying[i] = len(fullsat)
        
        if (len(current_population_truesat) > 0):
            current_population_truesat_unique = np.unique(current_population_truesat, axis = 0)
            for k2 in range(len(current_population_truesat_unique)):
                fullsat.append(current_population_truesat_unique[k2])
                           
        if (i != 0): #NEW
            previous_pareto_front = pareto_front_dict[nfe_array[i-1]].tolist() #NEW
            for k in range(len(previous_pareto_front)):
                current_population_unique.append(previous_pareto_front[k])
                
            #previous_pareto_constr_aggr = []
            #previous_pareto_truesat = []
            for design in previous_pareto_front:
                des_index = pen_obj1_constr_sorted.index(design[0])
                des_constr = []
                des_constr.append(1 - get_array_element(feas_scores_sorted, des_index))
                des_constr.append(1 - get_array_element(conn_scores_sorted, des_index))
                if (not artery_problem):
                     des_constr.append(get_array_element(stiffrat_vals_sorted, des_index))
                 
                constr_aggr = np.sum(np.abs(des_constr))
                current_population_constr_aggr.append(constr_aggr)
                # if (constr_aggr == 0):
                #      true_obj1, true_obj2 = get_true_objectives(true_obj1_sorted, true_obj2_sorted, des_index)
                #      current_population_truesat.append([-true_obj1, true_obj2])
                     
                previous_pareto_front_truesat = pareto_front_truesat_dict[nfe_array[i-1]].tolist()
                for k in range(len(previous_pareto_front_truesat)):
                    if (previous_pareto_front_truesat[k][0] == -1e4):
                        continue
                    else:
                        current_population_truesat.append(previous_pareto_front_truesat[k])
                     
        
        
        ######current_population_unique = np.unique(current_population, axis=0)
        
        #if ("AOS - Orient\\" in csv_filepath) and ("emoea_16" in csv_filepath):
            #print("stop")
        
        # current_population_constr_aggr = []
        # current_population_truesat = []
        # for design in current_population_unique:
        #     des_index = pen_obj1_constr_sorted.index(design[0])
        #     des_constr = []
        #     des_constr.append(1 - get_array_element(feas_scores_sorted, des_index))
        #     des_constr.append(1 - get_array_element(conn_scores_sorted, des_index))
        #     if (not artery_problem):
        #         des_constr.append(get_array_element(stiffrat_vals_sorted, des_index))
            
        #     constr_aggr = np.sum(np.abs(des_constr))
        #     current_population_constr_aggr.append(constr_aggr)
        #     if (constr_aggr == 0):
        #         true_obj1, true_obj2 = get_true_objectives(true_obj1_sorted, true_obj2_sorted, des_index)
        #         current_population_truesat.append([-true_obj1, true_obj2])
                
        # if (i != 0): #NEW
        #     previous_pareto_front_truesat = pareto_front_truesat_dict[nfe_array[i-1]].tolist()
        #     for k in range(len(previous_pareto_front_truesat)):
        #         if (previous_pareto_front_truesat[k][0] == -1e4):
        #             continue
        #         else:
        #             current_population_truesat.append(previous_pareto_front_truesat[k])
                
        
        if (not current_population_truesat):
            # purposefully add bad values so that compute_hv function can be used
            current_population_truesat.append([-1e4, 1.0])
            current_population_truesat_unique = np.unique(current_population_truesat, axis = 0)
            #num_fully_satisfying[i] = 0
            current_population_sat_constr_aggr = [5]
        else:
            current_population_truesat_unique = np.unique(current_population_truesat, axis = 0)
            #num_fully_satisfying[i] = len(current_population_truesat_unique)
            current_population_sat_constr_aggr = list(np.zeros(len(current_population_truesat_unique)))
        
        if (not fullsat):
            num_fully_satisfying[i] = 0
        else:
            fullsat_unique = np.unique(fullsat, axis = 0).tolist()
            num_fully_satisfying[i] = len(fullsat_unique)
            
        #if (num_fully_satisfying[i] < num_fully_satisfying[i-1]) and (i != 0):
            #print("exception")
        
        current_pareto_front_all = compute_pareto_front(current_population_unique, current_population_constr_aggr)
        #current_pareto_front = list(set(current_pareto_front_all))
        current_pareto_front = np.unique(current_pareto_front_all, axis=0)
        
        current_pareto_front_truesat_all = compute_pareto_front(current_population_truesat_unique, current_population_sat_constr_aggr)
        current_pareto_front_truesat = np.unique(current_pareto_front_truesat_all, axis=0)
    
        current_pareto_feas_scores = []
        current_pareto_conn_scores = []
        if (not artery_problem):
            current_pareto_stiffrat_vals = []
        #current_pareto_partcoll_scores = []
        #current_pareto_nodprop_scores = []
        #current_pareto_orient_scores = []
        #current_pareto_designs = []
        current_pareto_true_obj = []
        current_pareto_truesat_obj = []
        for pareto_design in current_pareto_front:
            design_index = pen_obj1_constr_sorted.index(pareto_design[0])
            design_feas_score = get_array_element(feas_scores_sorted, design_index)
            design_conn_score = get_array_element(conn_scores_sorted, design_index)
            if (not artery_problem):
                design_stiffrat_val = get_array_element(stiffrat_vals_sorted, design_index)
            #design_partcoll_score = get_array_element(partcoll_scores_sorted, design_index)
            #design_nodprop_score = get_array_element(nodprop_scores_sorted, design_index)
            #design_orient_score = get_array_element(orient_scores_sorted, design_index)
            current_pareto_feas_scores.append(design_feas_score)
            current_pareto_conn_scores.append(design_conn_score)
            if (not artery_problem):
                current_pareto_stiffrat_vals.append(design_stiffrat_val)
            #current_pareto_partcoll_scores.append(design_partcoll_score)
            #current_pareto_nodprop_scores.append(design_nodprop_score)
            #current_pareto_orient_scores.append(design_orient_score)
            #current_pareto_designs.append(get_array_element(designs_sorted, design_index))
            
            #true_obj1, true_obj2 = compute_true_objectives(pareto_design[0], pareto_design[1], design_feas_score, design_stab_score, fib_stif, intpen_feas_bool, intpen_stab_bool)
            true_obj1, true_obj2 = get_true_objectives(true_obj1_sorted, true_obj2_sorted, design_index)
            
            current_pareto_true_obj.append([true_obj1, true_obj2])
        
        pareto_front_dict[nfe_val] = current_pareto_front
        #pareto_front_feas_dict[nfe_val] = current_pareto_feas_scores
        #pareto_front_conn_dict[nfe_val] = current_pareto_conn_scores
        #pareto_front_stiffrat_dict[nfe_val] = current_pareto_stiffrat_vals
        #pareto_front_partcoll_dict[nfe_val] = current_pareto_partcoll_scores
        #pareto_front_nodprop_dict[nfe_val] = current_pareto_nodprop_scores
        #pareto_front_orient_dict[nfe_val] = current_pareto_orient_scores
        #pareto_front_designs_dict[nfe_val] = current_pareto_designs
        pareto_front_true_dict[nfe_val] = current_pareto_true_obj
        pareto_front_truesat_dict[nfe_val] = current_pareto_front_truesat
        
        pf_nfeval_obj1 = [row[0] for row in current_pareto_front]
        pf_nfeval_obj2 = [row[1] for row in current_pareto_front]
        
        pf_true_nfeval_obj1 = [row[0] for row in current_pareto_true_obj]
        pf_true_nfeval_obj2 = [row[1] for row in current_pareto_true_obj]
        
        pf_truesat_nfeval_obj1 = [row[0] for row in current_pareto_front_truesat]
        pf_truesat_nfeval_obj2 = [row[1] for row in current_pareto_front_truesat]
        
        pf_normalize_max_obj1.append(np.max(pf_nfeval_obj1))
        pf_normalize_min_obj1.append(np.min(pf_nfeval_obj1))
        pf_normalize_max_obj2.append(np.max(pf_nfeval_obj2))
        pf_normalize_min_obj2.append(np.min(pf_nfeval_obj2))
        
        pf_true_normalize_max_obj1.append(np.max(pf_true_nfeval_obj1))
        pf_true_normalize_min_obj1.append(np.min(pf_true_nfeval_obj1))
        pf_true_normalize_max_obj2.append(np.max(pf_true_nfeval_obj2))
        pf_true_normalize_min_obj2.append(np.min(pf_true_nfeval_obj2))
        
        pf_truesat_normalize_max_obj1.append(np.max(pf_truesat_nfeval_obj1))
        pf_truesat_normalize_min_obj1.append(np.min(pf_truesat_nfeval_obj1))
        pf_truesat_normalize_max_obj2.append(np.max(pf_truesat_nfeval_obj2))
        pf_truesat_normalize_min_obj2.append(np.min(pf_truesat_nfeval_obj2))
    
        #nonzero_feas_scores = True in (feas_score > 0.1 for feas_score in current_pareto_feas_scores)
        #if (nonzero_feas_scores):
            #if (not nfe_jump_recorded):
                #jump_nfe = i
                #nfe_jump_recorded = True
        
            #pareto_obj1s = [pareto_design[0] for pareto_design in current_pareto_front]
            #pareto_obj2s = [pareto_design[1] for pareto_design in current_pareto_front]
        
            #if (np.max(pareto_obj1s) > obj1_normalize_max_afterjump):
                #obj1_normalize_max_afterjump = np.max(pareto_obj1s)
        
            #if (np.max(pareto_obj2s) > obj2_normalize_max_afterjump):
                #obj2_normalize_max_afterjump = np.max(pareto_obj2s)
        
            #if (np.min(pareto_obj1s) < obj1_normalize_min_afterjump):
                #obj1_normalize_min_afterjump = np.min(pareto_obj1s)
        
            #if (np.min(pareto_obj2s) < obj2_normalize_min_afterjump):
                #obj2_normalize_min_afterjump = np.min(pareto_obj2s)

    ### Computing obj_normalize_fullrun, obj_normalize_afterjump and obj_normalized_true_fullrun using the entire run 
    #obj_normalize_max_fullrun = [np.max(pen_obj1_constr_sorted), np.max(pen_obj2_constr_sorted)]
    #obj_normalize_min_fullrun = [np.min(pen_obj1_constr_sorted), np.min(pen_obj2_constr_sorted)]

    #obj_true_normalize_max_fullrun = [np.max(true_obj1_sorted), np.max(true_obj2_sorted)]
    #obj_true_normalize_min_fullrun = [np.min(true_obj1_sorted), np.min(true_obj2_sorted)]

    # shift to above for loop if afterjump condition is used
    #obj1_normalize_max_afterjump = 0
    #obj1_normalize_min_afterjump = 0
    #obj2_normalize_max_afterjump = 0
    #obj2_normalize_min_afterjump = 0
    
    #obj_normalize_max_afterjump = [obj1_normalize_max_afterjump, obj2_normalize_max_afterjump]
    #obj_normalize_min_afterjump = [obj1_normalize_min_afterjump, obj2_normalize_min_afterjump]
    
    #obj_normalize_fullrun = [obj_normalize_min_fullrun, obj_normalize_max_fullrun]
    #obj_normalize_afterjump = [obj_normalize_min_afterjump, obj_normalize_max_afterjump]
    #obj_normalize_true_fullrun = [obj_true_normalize_min_fullrun, obj_true_normalize_max_fullrun]
    
    ### Computing obj_normalize_fullrun, obj_normalize_afterjump and obj_normalized_true_fullrun using the pareto fronts    
    obj_normalize_max_fullrun = [np.max(pf_normalize_max_obj1), np.max(pf_normalize_max_obj2)]
    obj_normalize_min_fullrun = [np.min(pf_normalize_min_obj1), np.min(pf_normalize_min_obj2)]
    
    obj_true_normalize_max_fullrun = [np.max(pf_true_normalize_max_obj1), np.max(pf_true_normalize_max_obj2)]
    obj_true_normalize_min_fullrun = [np.min(pf_true_normalize_min_obj1), np.min(pf_true_normalize_min_obj2)]
    
    obj_truesat_normalize_max_fullrun = [np.max(pf_truesat_normalize_max_obj1), np.max(pf_truesat_normalize_max_obj2)]
    obj_truesat_normalize_min_fullrun = [np.min(pf_truesat_normalize_min_obj1), np.min(pf_truesat_normalize_min_obj2)]

    obj_normalize_fullrun = [obj_normalize_min_fullrun, obj_normalize_max_fullrun]
    obj_normalize_true_fullrun = [obj_true_normalize_min_fullrun, obj_true_normalize_max_fullrun]
    obj_normalize_truesat_fullrun = [obj_truesat_normalize_min_fullrun, obj_truesat_normalize_max_fullrun]
    
    return num_fully_satisfying, pareto_front_dict, pareto_front_true_dict, pareto_front_truesat_dict, obj_normalize_fullrun, obj_normalize_true_fullrun, obj_normalize_truesat_fullrun, max_func_evals

#### Compute overall normalization objectives for single case study/all compared case studies (from the complete runs)
def compute_overall_norm_objs(objs_normalization_full, objs_normalization_true, objs_normalization_truesat):
    # Each input is a dictionary with key as the case study/run string and value as the corresponding 2D array
    
    obj1_max_full_allcases = np.zeros(len(objs_normalization_full))
    obj1_min_full_allcases = np.zeros(len(objs_normalization_full))
    obj2_max_full_allcases = np.zeros(len(objs_normalization_full))
    obj2_min_full_allcases = np.zeros(len(objs_normalization_full))
    obj1_max_true_allcases = np.zeros(len(objs_normalization_true))
    obj1_min_true_allcases = np.zeros(len(objs_normalization_true))
    obj2_max_true_allcases = np.zeros(len(objs_normalization_true))
    obj2_min_true_allcases = np.zeros(len(objs_normalization_true))
    obj1_max_truesat_allcases = np.zeros(len(objs_normalization_truesat))
    obj1_min_truesat_allcases = np.zeros(len(objs_normalization_truesat))
    obj2_max_truesat_allcases = np.zeros(len(objs_normalization_truesat))
    obj2_min_truesat_allcases = np.zeros(len(objs_normalization_truesat))
    
    i = 0
    for key in objs_normalization_full:
        current_objs_norm_full = objs_normalization_full[key]
        
        obj1_max_full_allcases[i] = current_objs_norm_full[1][0]
        obj2_max_full_allcases[i] = current_objs_norm_full[1][1]
        obj1_min_full_allcases[i] = current_objs_norm_full[0][0]
        obj2_min_full_allcases[i] = current_objs_norm_full[0][1]
        i += 1
        
    i = 0
    for key2 in objs_normalization_true:
        current_objs_norm_true = objs_normalization_true[key2]
        
        obj1_max_true_allcases[i] = current_objs_norm_true[1][0]
        obj2_max_true_allcases[i] = current_objs_norm_true[1][1]
        obj1_min_true_allcases[i] = current_objs_norm_true[0][0]
        obj2_min_true_allcases[i] = current_objs_norm_true[0][1]
        i += 1
        
    i = 0
    for key3 in objs_normalization_truesat:
        current_objs_norm_truesat = objs_normalization_truesat[key3]
        
        obj1_max_truesat_allcases[i] = current_objs_norm_truesat[1][0]
        obj2_max_truesat_allcases[i] = current_objs_norm_truesat[1][1]
        obj1_min_truesat_allcases[i] = current_objs_norm_truesat[0][0]
        obj2_min_truesat_allcases[i] = current_objs_norm_truesat[0][1]
        i += 1
        
    obj1_min_full_overall = np.min(obj1_min_full_allcases)
    obj2_min_full_overall = np.min(obj2_min_full_allcases)
    obj1_max_full_overall = np.max(obj1_max_full_allcases)
    obj2_max_full_overall = np.max(obj2_max_full_allcases)
    
    obj1_min_true_overall = np.min(obj1_min_true_allcases)
    obj2_min_true_overall = np.min(obj2_min_true_allcases)
    obj1_max_true_overall = np.max(obj1_max_true_allcases)
    obj2_max_true_overall = np.max(obj2_max_true_allcases)
    
    obj1_min_truesat_overall = np.min(obj1_min_truesat_allcases)
    obj2_min_truesat_overall = np.min(obj2_min_truesat_allcases)
    obj1_max_truesat_overall = np.max(obj1_max_truesat_allcases)
    obj2_max_truesat_overall = np.max(obj2_max_truesat_allcases)
            
    obj_norm_full_overall = [[obj1_min_full_overall, obj2_min_full_overall], [obj1_max_full_overall, obj2_max_full_overall]]
    obj_norm_true_overall = [[obj1_min_true_overall, obj2_min_true_overall], [obj1_max_true_overall, obj2_max_true_overall]]    
    obj_norm_truesat_overall = [[obj1_min_truesat_overall, obj2_min_truesat_overall], [obj1_max_truesat_overall, obj2_max_truesat_overall]]    
    
    return obj_norm_full_overall, obj_norm_true_overall, obj_norm_truesat_overall

#### Compute hypervolume arrays from copmuted pareto fronts and normalization constants
def compute_hv_arrays_from_csv_data(pf_dict, pf_true_dict, pf_truesat_dict, obj_norm_full, obj_norm_true_full, obj_norm_truesat_full, max_fun_evals):
    obj_norm_min_full = obj_norm_full[0]
    obj_norm_max_full = obj_norm_full[1]
    obj_norm_true_min_full = obj_norm_true_full[0]
    obj_norm_true_max_full = obj_norm_true_full[1]
    obj_norm_truesat_min_full = obj_norm_truesat_full[0]
    obj_norm_truesat_max_full = obj_norm_truesat_full[1]

    ## Normalize the pareto front objectives and compute the hypervolume
    hypervol_full_dict = []
    hypervol_true_full_dict = []
    hypervol_truesat_full_dict = []

    for nfe_val in nfe_array:
        #print('iter = ' + str(nfe_val))
    
        current_pareto_front = pf_dict[nfe_val]
        current_true_pareto_front = pf_true_dict[nfe_val]
        current_truesat_pareto_front = pf_truesat_dict[nfe_val]
        current_pf_normalized = []
        current_pf_true_normalized = []
        current_pf_truesat_normalized = []
        for pareto_design in current_pareto_front:
            obj1_normalized = (pareto_design[0] - obj_norm_min_full[0])/(obj_norm_max_full[0] - obj_norm_min_full[0])
            obj2_normalized = (pareto_design[1] - obj_norm_min_full[1])/(obj_norm_max_full[1] - obj_norm_min_full[1])
            current_pf_normalized.append([obj1_normalized, obj2_normalized])
            
        for pareto_design_true in current_true_pareto_front:
            obj1_true_normalized = (obj_norm_true_max_full[0] - pareto_design_true[0])/(obj_norm_true_max_full[0] - obj_norm_true_min_full[0])
            obj2_true_normalized = (pareto_design_true[1] - obj_norm_true_min_full[1])/(obj_norm_true_max_full[1] - obj_norm_true_min_full[1])
            current_pf_true_normalized.append([obj1_true_normalized, obj2_true_normalized])
            
        for pareto_design_truesat in current_truesat_pareto_front:
            obj1_truesat_normalized = (pareto_design_truesat[0] - obj_norm_truesat_min_full[0])/(obj_norm_truesat_max_full[0] - obj_norm_truesat_min_full[0])
            obj2_truesat_normalized = (pareto_design_truesat[1] - obj_norm_truesat_min_full[1])/(obj_norm_truesat_max_full[1] - obj_norm_truesat_min_full[1])
            current_pf_truesat_normalized.append([obj1_truesat_normalized, obj2_truesat_normalized])
            
        current_hv = compute_hv(current_pf_normalized)
        hypervol_full_dict.append([nfe_val, current_hv])
        
        current_hv_true = compute_hv(current_pf_true_normalized)
        hypervol_true_full_dict.append([nfe_val, current_hv_true])
        
        current_hv_truesat = compute_hv(current_pf_truesat_normalized)
        hypervol_truesat_full_dict.append([nfe_val, current_hv_truesat])
        
    return hypervol_full_dict, hypervol_true_full_dict, hypervol_truesat_full_dict

### Compute array of NFE values for reaching threshold hypervolume for different runs of a particular case
def compute_nfe_hypervolume_attained(hv_dict, hv_threshold):
    #hv_threshold = Threshold HV value to reach, user parameter
    n_runs = len(hv_dict)
    nfe_hv_attained = []
    for key in hv_dict:
        hv_array_run = hv_dict[key]
        nfe_array_run = [hv_array[0] for hv_array in hv_array_run]
        hv_val_array = [hv_array[1] for hv_array in hv_array_run]
        nfe_hv_attained_run = nfe_array_run[-1] + 100
        for i in range(len(hv_val_array)):
            if (hv_val_array[i] >= hv_threshold):
                nfe_hv_attained_run = nfe_array_run[i]
                break
                
        #hv_val_diff = [np.abs(x - hv_threshold) for x in hv_val_array]
        #index = np.argmin(hv_val_diff)
        nfe_hv_attained.append(nfe_hv_attained_run)
        
    return nfe_hv_attained
    
### Plot fraction of runs attaining threshold hypervolume vs NFE
def plot_fraction_hypervolume_attained(nfe_hv_attained_dict, nfe_array, colour_array, casename_array, savefig_name, hv_threshold):
    fig1 = plt.figure()
    n_cases = len(nfe_hv_attained_dict)
    case_idx = 0
    for case_key in nfe_hv_attained_dict:
        nfe_hv_attained_case = nfe_hv_attained_dict[case_key]
        n_runs = len(nfe_hv_attained_case)
        frac_runs_hv_attained = np.zeros(len(nfe_array))
        for i in range(len(nfe_array)):
            idx_runs_hv_attained = [idx for idx, val in enumerate(nfe_hv_attained_case) if val <= nfe_array[i]]
            frac_runs_hv_attained[i] = len(idx_runs_hv_attained)/n_runs
            
        plt.plot(nfe_array, frac_runs_hv_attained, '-', color=colour_array[case_idx], label=casename_array[case_idx])
        case_idx += 1
    
    plt.xlabel('Number of Function Evaluations',fontsize=14)
    plt.ylabel('Fraction of runs HV $\geq$ ' + str(hv_threshold),fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.15), ncol=3, borderaxespad=0 ,prop={"size":12})
    plt.show()
    #fig1.savefig('frac_hv_attained_' + savefig_name + '.png', format='png')
    
def compute_hypervolume_stats(hypervols_dict):
    hv_dict_keys = list(hypervols_dict.keys())
    hv_dict_0 = hypervols_dict[hv_dict_keys[0]]
    nfe_array_0 = [hv_array[0] for hv_array in hv_dict_0]
    n_datapoints = len(nfe_array_0)
    hypervol_median = np.zeros(n_datapoints)
    hypervol_1q = np.zeros(n_datapoints)
    hypervol_3q = np.zeros(n_datapoints)
    for i in range(n_datapoints):
        hypervol_vals = []
        for key in hypervols_dict:
            hv_dict_j = hypervols_dict[key]
            hv_current_array = [hv_array[1] for hv_array in hv_dict_j]
            hypervol_vals.append(hv_current_array[i])
        hypervol_median[i] = statistics.median(hypervol_vals)
        #hypervol_median[i] = statistics.mean(hypervol_vals)
        hypervol_1q[i] = np.percentile(hypervol_vals, 25)
        #hypervol_1q[i] = hypervol_median[i] - statistics.stdev(hypervol_vals)
        hypervol_3q[i] = np.percentile(hypervol_vals, 75)
        #hypervol_3q[i] = hypervol_median[i] + statistics.stdev(hypervol_vals)
        
    return hypervol_median, hypervol_1q, hypervol_3q, nfe_array_0

def compute_mann_whitney_Uvals(artery_prob, hv_dict_allcases_allruns, nfe_array): # Wilcoxon Rank Sum Test
    ## hv_med_array_allcases is a dictionary of length = number of cases
    
    if artery_prob:
        ## For artery problem
        nfe_samples_array = [0, 250, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
        n_samples = len(nfe_samples_array)
    else:
        ## For truss problem
        n_samples = 5
        linspace_samples_array = np.linspace(0,1,n_samples)
        nfe_samples_array = np.floor(np.multiply(linspace_samples_array, nfe_array[-1]))
    
    nfe_samples_indices_array = np.zeros(n_samples)
    for i in range(len(nfe_samples_array)):
        nfe_samples_indices_array[i] = find_closest_index(nfe_samples_array[i], nfe_array)
        
    hv_samples_allcases_allruns = {}
    for j in range(len(hv_dict_allcases_allruns)):
        hv_dict_allruns_currentcase = hv_dict_allcases_allruns['case'+str(j)]
        hv_samples_allruns = {}
        for k in range(len(nfe_samples_indices_array)):
            hv_samples_nfe_allruns = np.zeros(len(hv_dict_allruns_currentcase))
            for m in range(len(hv_dict_allruns_currentcase)):
                hv_run = hv_dict_allruns_currentcase['run'+str(m)]
                hv_samples_nfe_allruns[m] = hv_run[int(nfe_samples_indices_array[k])][1]
            hv_samples_allruns['nfe:'+str(int(nfe_samples_array[k]))] = hv_samples_nfe_allruns
        hv_samples_allcases_allruns['case'+str(j)] = hv_samples_allruns
    
    cases_array = np.arange(len(hv_dict_allcases_allruns))
    case_combinations = list(combinations(cases_array,2))
    
    U_test_cases = {}
    
    for n in range(len(case_combinations)):
        case_string_x = 'case' + str(case_combinations[n][0])
        case_string_y = 'case' + str(case_combinations[n][1])
        
        hv_allruns_casex = hv_samples_allcases_allruns[case_string_x]
        hv_allruns_casey = hv_samples_allcases_allruns[case_string_y]
        
        U_test_cases_allnfes = {}
        for p in range(len(nfe_samples_indices_array)):
            hv_samples_nfe_casex = hv_allruns_casex['nfe:'+str(int(nfe_samples_array[p]))]
            hv_samples_nfe_casey = hv_allruns_casey['nfe:'+str(int(nfe_samples_array[p]))]
            
            U1, p_val = mannwhitneyu(hv_samples_nfe_casex, hv_samples_nfe_casey, alternative='less')
            t_val, p_val_t = ttest_ind(hv_samples_nfe_casex, hv_samples_nfe_casey, equal_var=False, alternative='less')
            
            U2 = len(hv_samples_nfe_casex)*len(hv_samples_nfe_casey) - U1
            
            U_test = np.min(np.array([U1, U2]))
            
            U_test_cases_allnfes['nfe:'+str(int(nfe_samples_array[p]))] = [U_test, p_val, p_val_t]
        
        dict_key = case_string_x + ' and ' + case_string_y
        U_test_cases[dict_key] = U_test_cases_allnfes
        
    return U_test_cases

def plot_number_fully_satisfying_designs(num_fullsat_allcases, nfe_array, colour_array, alpha_array, casename_array):
    fig = plt.figure()
    ## Compute stats of number of fully satisfying designs 
    num_sat_allcases_keys = list(num_fullsat_allcases.keys())
    num_cases = len(num_fullsat_allcases)
    for i in range(num_cases):
        key = num_sat_allcases_keys[i]
        num_sat_allruns = num_fullsat_allcases[key]
        run_keys = list(num_sat_allruns.keys())
        num_sat_mean = np.zeros(len(num_sat_allruns[run_keys[0]]))
        num_sat_upstd = np.zeros(len(num_sat_allruns[run_keys[0]]))
        num_sat_downstd = np.zeros(len(num_sat_allruns[run_keys[0]]))
        for j in range(len(num_sat_allruns[run_keys[0]])):
            num_sat_nfe_allruns = []
            for run_key in run_keys:
                num_sat_run = num_sat_allruns[run_key]
                num_sat_nfe_allruns.append(num_sat_run[j])
            #num_sat_mean[j] = statistics.mean(num_sat_nfe_allruns)
            num_sat_mean[j] = statistics.median(num_sat_nfe_allruns)
            #num_sat_upstd[j] = num_sat_mean[j] + 2*statistics.stdev(num_sat_nfe_allruns)
            num_sat_upstd[j] = np.percentile(num_sat_nfe_allruns, 75)
            #num_sat_downstd[j] = num_sat_mean[j] - 2*statistics.stdev(num_sat_nfe_allruns)        
            num_sat_downstd[j] = np.percentile(num_sat_nfe_allruns, 25)
        plt.plot(nfe_array, num_sat_mean, '-', color=colour_array[i], label=casename_array[i])
        plt.fill_between(nfe_array, num_sat_downstd, num_sat_upstd, color=colour_array[i], alpha=alpha_array[i])
    plt.xlabel('Number of Function Evaluations',fontsize=14)
    plt.ylabel('Number of Fully Satisfying Designs',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.title(plot_title)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.15), ncol=3, borderaxespad=0, prop={"size":12})
    plt.show()

def plot_hypervolume_stats(hv_median_case, hv_1q_case, hv_3q_case, nfe_array, savefig_name):
    fig1 = plt.figure()
    plt.plot(nfe_array,hv_median_case, 'b-', label='Median')
    plt.plot(nfe_array,hv_1q_case, 'r-', label='1st Quartile')
    plt.plot(nfe_array,hv_3q_case, 'g-', label='3rd Quartile')
    plt.xlabel('Number of Function Evaluations')
    plt.ylabel('Hypervolume')
    plt.title('Averaged Hypervolume vs NFE')
    plt.legend(loc='lower right')
    plt.show()
    #fig1.savefig('HV_plot_averaged_' + savefig_name + '.png')
    
def plot_hypervolume_stats_case_allruns(hv_dict_allruns, nfe_array, case_color, casename):
    fig = plt.figure()
    number_runs = len(hv_dict_allruns)
    for i in range(number_runs):
        hv_vals = [x[1] for x in hv_dict_allruns['run'+str(i)]]
        plt.plot(nfe_array, hv_vals, '-', color=case_color)
    plt.xlabel('Number of Function Evaluations',fontsize=14)
    plt.ylabel('Hypervolume',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title('Hypervolume (only satisfying designs) - all runs for ' + casename)
    plt.show()
    
def plot_hypervolume_stats_allcases(hv_median_dict, hv_1q_dict, hv_3q_dict, nfe_array, colour_array, alpha_array, casename_array, plot_title, savefig_name):
    fig1 = plt.figure()
    number_cases = len(hv_median_dict)
    #print('n_cases')
    #print(number_cases)
    for i in range(number_cases):
        #print(print(marker_array[i]+'*'))
        plt.plot(nfe_array, hv_median_dict['case'+str(i)], '-', color=colour_array[i], label=casename_array[i])
        plt.fill_between(nfe_array, hv_1q_dict['case'+str(i)], hv_3q_dict['case'+str(i)], color=colour_array[i], alpha=alpha_array[i])
        
        #plt.plot(nfe_array, hv_1q_dict['case'+str(i)], '--', color=colour_array[i])#, label=casename_array[i]+' 1st Quartile')
        #plt.plot(nfe_array, hv_3q_dict['case'+str(i)], '--', color=colour_array[i])#, label=casename_array[i]+' 3rd Quartile')
    plt.xlabel('Number of Function Evaluations',fontsize=14)
    plt.ylabel('Hypervolume',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.title(plot_title)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.15), ncol=3, borderaxespad=0, prop={"size":12})
    plt.show()
    #fig1.savefig('HV_plot_averaged_' + savefig_name + '.png', format='png')
    
#### Define functions to compute and plot hypervolume for single case and all cases
def hypervolume_computation_single_case(choice_model, credit_assign, case_booleans, prob_artery, sidenum, run_nums, case_name):
    ## Computing the pareto fronts and normalization objectives for each run
    obj_norm_allruns = {}
    obj_norm_true_allruns = {}
    obj_norm_truesat_allruns = {}
    pf_allruns = {}
    pf_true_allruns = {}
    pf_truesat_allruns = {}
    num_sat_allruns = {}
    max_f_evals_allruns = np.zeros(run_nums)
    for i in range(run_nums):
        print('Computing Pareto Fronts for run ' + str(i))
        current_csvpath = get_csv_filepath_material(case_booleans[:4], case_booleans[4:8], case_booleans[8:12], case_booleans[12:16], prob_artery, choice_model, credit_assign, i)
        heur_intpen_constr = [case_booleans[0], case_booleans[4], case_booleans[8], case_booleans[12]]
        num_sat_run, pf_dict_i, pf_true_dict_i, pf_truesat_dict_i, obj_norm_full_i, obj_norm_true_i, obj_norm_truesat_i, max_fun_evals_i = extract_data_from_csv(current_csvpath, prob_artery, heur_intpen_constr, sidenum)
        pf_allruns['run'+str(i)] = pf_dict_i
        pf_true_allruns['run'+str(i)] = pf_true_dict_i
        pf_truesat_allruns['run'+str(i)] = pf_truesat_dict_i
        obj_norm_allruns['run'+str(i)] = obj_norm_full_i
        obj_norm_true_allruns['run'+str(i)] = obj_norm_true_i
        obj_norm_truesat_allruns['run'+str(i)] = obj_norm_truesat_i
        num_sat_allruns['run'+str(i)] = num_sat_run
        max_f_evals_allruns[i] = max_fun_evals_i
    
    #print('pf_allruns')
    #print(pf_allruns)
    #print('\n')
    #print('pf_true_allruns')
    #print(pf_true_allruns)
    #print('\n')
    ## Use computed normalization objectives and find the overall normalization objectives across all runs
    print('Computing overall normalization constants')
    norm_objs_full_overall, norm_objs_true_overall, norm_objs_truesat_overall = compute_overall_norm_objs(obj_norm_allruns, obj_norm_true_allruns, obj_norm_truesat_allruns)
    #print('norm_objs_full_overall')
    #print(norm_objs_full_overall)
    #print('\n')
    #print('norm_objs_true_overall')
    #print(norm_objs_true_overall)
    #print('\n')
    
    ## Compute Hypervolume values for each run
    hv_dict_allruns = {}
    #hv_dict_aj_allruns = {}
    hv_dict_true_allruns = {}
    hv_dict_truesat_allruns = {}
    for j in range(run_nums):
        print('Computing hypervolume values for run ' + str(j))
        hv_dict_j, hv_dict_true_j, hv_dict_truesat_j = compute_hv_arrays_from_csv_data(pf_allruns['run'+str(j)], pf_true_allruns['run'+str(j)], pf_truesat_allruns['run'+str(j)], norm_objs_full_overall, norm_objs_true_overall, norm_objs_truesat_overall, max_f_evals_allruns[j])
        hv_dict_allruns['run'+str(j)] = hv_dict_j
        hv_dict_true_allruns['run'+str(j)] = hv_dict_true_j
        hv_dict_truesat_allruns['run'+str(j)] = hv_dict_truesat_j
        
    #print('hv_dict_allruns')
    #print(hv_dict_allruns)
    #print('hv_dict_true_allruns')
    #print(hv_dict_true_allruns)
        
    ## Plotting
    print('Plotting')
    hv_median_all, hv_1q_all, hv_3q_all, nfe_array = compute_hypervolume_stats(hv_dict_allruns)
    plot_hypervolume_stats(hv_median_all, hv_1q_all, hv_3q_all, nfe_array, case_name+'_full')

    ## Plot HVs for true objectives
    hv_true_median_all, hv_true_1q_all, hv_true_3q_all, nfe_array_true = compute_hypervolume_stats(hv_dict_true_allruns)
    plot_hypervolume_stats(hv_true_median_all, hv_true_1q_all, hv_true_3q_all, nfe_array_true, case_name+'_true')
    
    ## Plot HVs for truesat objectives
    hv_truesat_median_all, hv_truesat_1q_all, hv_truesat_3q_all, nfe_array_truesat = compute_hypervolume_stats(hv_dict_truesat_allruns)
    plot_hypervolume_stats(hv_truesat_median_all, hv_truesat_1q_all, hv_truesat_3q_all, nfe_array_truesat, case_name+'_truesat')
    
def hypervolume_computation_all_cases(choice_model, credit_assign, case_bools_dict, prob_artery, sidenum, run_nums, marker_colours, alpha_vals, case_names, hv_thresh):
    num_cases = len(case_bools_dict) # number of cases to compare 

    ## Computing the pareto fronts and normalization objectives for each run in each case
    pf_allcases = {}
    pf_true_allcases = {}
    pf_truesat_allcases = {}
    obj_norm_allcasesandruns = {}
    obj_norm_true_allcasesandruns = {}
    obj_norm_truesat_allcasesandruns = {}
    num_sat_allcases = {}
    max_f_evals_allcases = {}
    for i in range(num_cases):
        print('Computing Pareto Fronts for runs in Case '+str(i))
        current_case_bools = case_bools_dict['case'+str(i+1)]
        pf_allruns_i = {}
        pf_true_allruns_i = {}
        pf_truesat_allruns_i = {}
        num_sat_allruns = {}
        max_f_evals_allruns = np.zeros(run_nums)
        for j in range(run_nums):
            print('Run '+str(j))
            current_csvpath = get_csv_filepath_material(current_case_bools[:4], current_case_bools[4:8], current_case_bools[8:12], current_case_bools[12:16], prob_artery, choice_model, credit_assign, j)
            heur_intpen_constr = [current_case_bools[0], current_case_bools[4], current_case_bools[8], current_case_bools[12]]
            num_sat_run, pf_dict_j, pf_true_dict_j, pf_truesat_dict_j, obj_norm_full_j, obj_norm_true_j, obj_norm_truesat_j, max_fun_evals_j = extract_data_from_csv(current_csvpath, prob_artery, heur_intpen_constr, sidenum)
            pf_allruns_i['run'+str(j)] = pf_dict_j
            pf_true_allruns_i['run'+str(j)] = pf_true_dict_j
            pf_truesat_allruns_i['run'+str(j)] = pf_truesat_dict_j
            obj_norm_allcasesandruns['case'+str(i+1)+'run'+str(j)] = obj_norm_full_j
            obj_norm_true_allcasesandruns['case'+str(i+1)+'run'+str(j)] = obj_norm_true_j
            obj_norm_truesat_allcasesandruns['case'+str(i+1)+'run'+str(j)] = obj_norm_truesat_j
            num_sat_allruns['run'+str(j)] = num_sat_run
            max_f_evals_allruns[j] = max_fun_evals_j
        pf_allcases['case'+str(i+1)] = pf_allruns_i
        pf_true_allcases['case'+str(i+1)] = pf_true_allruns_i
        pf_truesat_allcases['case'+str(i+1)] = pf_truesat_allruns_i
        num_sat_allcases['case'+str(i+1)] = num_sat_allruns
        max_f_evals_allcases['case'+str(i+1)] = max_f_evals_allruns
    
    ## Use computed normalization objectives and find the overall normalization objectives across all runs and cases
    print('Computing overall normalization constants')
    norm_objs_full_overall, norm_objs_true_overall, norm_objs_truesat_overall = compute_overall_norm_objs(obj_norm_allcasesandruns, obj_norm_true_allcasesandruns, obj_norm_truesat_allcasesandruns)
    
    ## Compute Hypervolume values for each run in each case
    hv_dict_allcases = {}
    hv_dict_true_allcases = {}
    hv_dict_truesat_allcases = {}
    hv_dict_median_allcases = {}
    hv_dict_1q_allcases = {}
    hv_dict_3q_allcases = {}
    #hv_dict_aj_median_allcases = {}
    #hv_dict_aj_1q_allcases = {}
    #hv_dict_aj_3q_allcases = {}
    hv_dict_true_median_allcases = {}
    hv_dict_true_1q_allcases = {}
    hv_dict_true_3q_allcases = {}
    hv_dict_truesat_median_allcases = {}
    hv_dict_truesat_1q_allcases = {}
    hv_dict_truesat_3q_allcases = {}
    nfe_array_hv_attained_dict = {}
    for i in range(num_cases):
        print('Computing hypervolume values for runs in Case '+str(i))
        pfs_case_i = pf_allcases['case'+str(i+1)]
        pfs_true_case_i = pf_true_allcases['case'+str(i+1)]
        pfs_truesat_case_i = pf_truesat_allcases['case'+str(i+1)]
        max_func_evals_i = max_f_evals_allcases['case'+str(i+1)]
        hv_dict_allruns = {}
        hv_dict_true_allruns = {}
        hv_dict_truesat_allruns = {}
        for j in range(run_nums):
            print('Run '+str(j))
            hv_dict_j, hv_dict_true_j, hv_dict_truesat_j = compute_hv_arrays_from_csv_data(pfs_case_i['run'+str(j)], pfs_true_case_i['run'+str(j)], pfs_truesat_case_i['run'+str(j)], norm_objs_full_overall, norm_objs_true_overall, norm_objs_truesat_overall, max_func_evals_i[j])
            hv_dict_allruns['run'+str(j)] = hv_dict_j
            hv_dict_true_allruns['run'+str(j)] = hv_dict_true_j
            hv_dict_truesat_allruns['run'+str(j)] = hv_dict_truesat_j
            
        hv_dict_allcases['case'+str(i)] = hv_dict_allruns
        hv_dict_true_allcases['case'+str(i)] = hv_dict_true_allruns
        hv_dict_truesat_allcases['case'+str(i)] = hv_dict_truesat_allruns
        
        print('Plotting hypervolume plots for all runs of case ' + str(i))
        plot_hypervolume_stats_case_allruns(hv_dict_truesat_allruns, nfe_array, marker_colours[i], case_names[i])
        #plot_hypervolume_stats_case_allruns(hv_dict_allruns, nfe_array, case_names[i])
        
        print('Computing array of NFE for attaining threshold hypervolume')
        nfe_hv_attained_case = compute_nfe_hypervolume_attained(hv_dict_truesat_allruns, hv_thresh)
        nfe_array_hv_attained_dict['case'+str(i)] = nfe_hv_attained_case
                
        print('Computing hypervolume stats')
        hv_med_i, hv_1q_i, hv_3q_i, nfe_array_i = compute_hypervolume_stats(hv_dict_allruns)
        hv_med_true_i, hv_1q_true_i, hv_3q_true_i, nfe_array_true_i = compute_hypervolume_stats(hv_dict_true_allruns)
        hv_med_truesat_i, hv_1q_truesat_i, hv_3q_truesat_i, nfe_array_truesat_i = compute_hypervolume_stats(hv_dict_truesat_allruns)
        
        hv_dict_median_allcases['case'+str(i)] = hv_med_i
        hv_dict_1q_allcases['case'+str(i)] = hv_1q_i
        hv_dict_3q_allcases['case'+str(i)] = hv_3q_i
        hv_dict_true_median_allcases['case'+str(i)] = hv_med_true_i
        hv_dict_true_1q_allcases['case'+str(i)] = hv_1q_true_i
        hv_dict_true_3q_allcases['case'+str(i)] = hv_3q_true_i
        hv_dict_truesat_median_allcases['case'+str(i)] = hv_med_truesat_i
        hv_dict_truesat_1q_allcases['case'+str(i)] = hv_1q_truesat_i
        hv_dict_truesat_3q_allcases['case'+str(i)] = hv_3q_truesat_i
        
    print('Computing Wilcoxon Test Statistics')
    U_test_dict = compute_mann_whitney_Uvals(prob_artery, hv_dict_allcases, nfe_array_i)
    U_test_truesat_dict = compute_mann_whitney_Uvals(prob_artery, hv_dict_truesat_allcases, nfe_array_truesat_i)
    
    #return nfe_array_hv_attained_dict, hv_dict_median_allcases, hv_dict_1q_allcases, hv_dict_3q_allcases, hv_dict_true_median_allcases, hv_dict_true_1q_allcases, hv_dict_true_3q_allcases, hv_dict_truesat_median_allcases, hv_dict_truesat_1q_allcases, hv_dict_truesat_3q_allcases, nfe_array_i
    return num_sat_allcases, nfe_array_hv_attained_dict, hv_dict_median_allcases, hv_dict_1q_allcases, hv_dict_3q_allcases, hv_dict_true_median_allcases, hv_dict_true_1q_allcases, hv_dict_true_3q_allcases, hv_dict_truesat_median_allcases, hv_dict_truesat_1q_allcases, hv_dict_truesat_3q_allcases, U_test_dict, U_test_truesat_dict, nfe_array_i
    
def plotting_all_cases(num_sat_allcases, nfe_hv_attained_dict, hv_dict_med_allcases, hv_dict_1stq_allcases, hv_dict_3rdq_allcases, hv_dict_true_med_allcases, hv_dict_true_1stq_allcases, hv_dict_true_3rdq_allcases, hv_dict_truesat_med_allcases, hv_dict_truesat_1stq_allcases, hv_dict_truesat_3rdq_allcases, nfe_array0, mark_colors, alphas, names_cases, hv_thresh):
    print('Plotting')
    plot_number_fully_satisfying_designs(num_sat_allcases, nfe_array0, mark_colors, alphas, names_cases)
    plot_fraction_hypervolume_attained(nfe_hv_attained_dict, nfe_array0, mark_colors, names_cases, 'allcases_full', hv_thresh)
    plot_hypervolume_stats_allcases(hv_dict_med_allcases, hv_dict_1stq_allcases, hv_dict_3rdq_allcases, nfe_array, mark_colors, alphas, names_cases, 'Hypervolume of Normalized Objectives', 'allcases_full')
    ## Plot HVs for hv_afterjump
    plot_hypervolume_stats_allcases(hv_dict_true_med_allcases, hv_dict_true_1stq_allcases, hv_dict_true_3rdq_allcases, nfe_array, mark_colors, alphas, names_cases, 'Hypervolume of True Objectives', 'allcases_true')
    plot_hypervolume_stats_allcases(hv_dict_truesat_med_allcases, hv_dict_truesat_1stq_allcases, hv_dict_truesat_3rdq_allcases, nfe_array, mark_colors, alphas, names_cases, 'Hypervolume of True Objectives for fully satisfying designs', 'allcases_truesat')
    
    
#######################################################################################################################################################################################################################################################################################################################################################################################################
    # PROGRAM OPERATION
#######################################################################################################################################################################################################################################################################################################################################################################################################




#### Comparing Simple E-MOEA with AOS - all heuristics and AOS - promising heuristics
model_used = 2 # 1 = Fibre stiffness, 2 = Truss stiffness, 3 = APDL Beam
sidenum = 3 # 3x3 node grid
cases_dict = {}
artery_problem = True
num_runs = 30 # number of runs for each case
threshold_hv = 0.65

credit_assignment = 2 # 0 -> offspring parent dominance, 1 -> set improvement dominance, 2 -> set contribution dominance

# bools = [int_pen_partcoll, AOS_partcoll, bias_init_partcoll, ACH_partcoll, int_pen_nodalprop, AOS_nodalprop, bias_init_nodalprop, ACH_nodalprop, int_pen_orient, AOS_orient, bias_init_orient, ACH_orient, int_pen_inters, AOS_inters, bias_init_inters, ACH_inters]
case1_bools = [False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False] # Simple E-MOEA
case2_bools = [False, True, False, False, False, True, False, False, False, True, False, False, False, True, False, False] #  AOS - all heuristics
if artery_problem:
    case3_bools = [False, False, False, False, False, False, False, False, False, True, False, False, False, True, False, False] #  AOS - Orientation and Intersection
else:
    case3_bools = [False, False, False, False, False, False, False, False, False, True, False, False, False, True, False, False] #  AOS - Orientation and Intersection

cases_dict['case1'] = case1_bools
cases_dict['case2'] = case2_bools
cases_dict['case3'] = case3_bools

line_colours = ['#000000','#E69F00','#56B4E9'] # black, yellow, blue 
casenames = ['Eps. MOEA','All heurs','Promising heurs']

alpha_values = [0.5,0.5,0.5] # change based on number of cases/visibility

#nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, nfe_array_1 = hypervolume_computation_all_cases(model_used, cases_dict, artery_problem, sidenum, num_runs, line_colours, alpha_values, casenames, threshold_hv)
num_fullysat_cases, nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, Uvals_test, Uvals_test_truesat, nfe_array_1 = hypervolume_computation_all_cases(model_used, credit_assignment, cases_dict, artery_problem, sidenum, num_runs, line_colours, alpha_values, casenames, threshold_hv)

## Mann Whitney U values
print("For optimization objectives")
print(Uvals_test)
print("\n")
print("For true objectives of fully satisfying designs")
print(Uvals_test_truesat)

casenames = ['Eps. MOEA.','All heurs','Promising heurs']
plotting_all_cases(num_fullysat_cases, nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, nfe_array_1, line_colours, alpha_values, casenames, threshold_hv)
#print(nfe_cdf_array)
#print(hv_dict_truesat_med_cases)



# #### Comparing Simple E-MOEA with Int Pen - PartColl, AOS - PartColl, Bias Init - PartColl and ACH - PartColl
# model_used = 2 # 1 = Fibre stiffness, 2 = Truss stiffness, 3 = APDL Beam
# sidenum = 3 # 3x3 node grid
# cases_dict = {}
# artery_problem = False
# num_runs = 30 # number of runs for each case
# threshold_hv = 0.58

# # bools = [int_pen_partcoll, AOS_partcoll, bias_init_partcoll, ACH_partcoll, int_pen_nodalprop, AOS_nodalprop, bias_init_nodalprop, ACH_nodalprop, int_pen_orient, AOS_orient, bias_init_orient, ACH_orient, int_pen_inters, AOS_inters, bias_init_inters, ACH_inters]
# case1_bools = [False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False] # Simple E-MOEA
# case2_bools = [True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False] # Simple E-MOEA
# case3_bools = [False, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False] #  
# case4_bools = [False, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False] # 
# case5_bools = [False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, False] # 

# cases_dict['case1'] = case1_bools
# cases_dict['case2'] = case2_bools
# cases_dict['case3'] = case3_bools
# cases_dict['case4'] = case4_bools
# cases_dict['case5'] = case5_bools

# line_colours = ['#000000','#E69F00','#56B4E9','#FF0E0E','#30FF0E'] # black, yellow, blue, red, green  
# casenames = ['Eps. MOEA','Int Pen - PartColl','AOS - PartColl','Bias Init - PartColl','ACH - PartColl']

# alpha_values = [0.5,0.5,0.5,0.5,0.5] # change based on number of cases/visibility

# #nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, nfe_array_1 = hypervolume_computation_all_cases(model_used, cases_dict, artery_problem, sidenum, num_runs, line_colours, alpha_values, casenames, threshold_hv)
# num_fullysat_cases, nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, Uvals_test, Uvals_test_truesat, nfe_array_1 = hypervolume_computation_all_cases(model_used, cases_dict, artery_problem, sidenum, num_runs, line_colours, alpha_values, casenames, threshold_hv)

# ## Mann Whitney U values
# print("For optimization objectives")
# print(Uvals_test)
# print("\n")
# print("For true objectives of fully satisfying designs")
# print(Uvals_test_truesat)

# casenames = ['Eps. MOEA','Int Pen - PartColl','AOS - PartColl','Bias Init - PartColl','ACH - PartColl']
# plotting_all_cases(num_fullysat_cases, nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, nfe_array_1, line_colours, alpha_values, casenames, threshold_hv)
# #print(nfe_cdf_array)
# #print(hv_dict_truesat_med_cases)



# #### Comparing Simple E-MOEA with Int Pen - NodalProp, AOS - NodalProp, Bias Init - NodalProp and ACH - NodalProp
# model_used = 2 # 1 = Fibre stiffness, 2 = Truss stiffness, 3 = APDL Beam
# sidenum = 3 # 3x3 node grid
# cases_dict = {}
# artery_problem = False
# num_runs = 30 # number of runs for each case
# threshold_hv = 0.58

# # bools = [int_pen_partcoll, AOS_partcoll, bias_init_partcoll, ACH_partcoll, int_pen_nodalprop, AOS_nodalprop, bias_init_nodalprop, ACH_nodalprop, int_pen_orient, AOS_orient, bias_init_orient, ACH_orient, int_pen_inters, AOS_inters, bias_init_inters, ACH_inters]
# case1_bools = [False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False] # Simple E-MOEA
# case2_bools = [False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False] # Simple E-MOEA
# case3_bools = [False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False] #  
# case4_bools = [False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False] # 
# case5_bools = [False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False] # 

# cases_dict['case1'] = case1_bools
# cases_dict['case2'] = case2_bools
# cases_dict['case3'] = case3_bools
# cases_dict['case4'] = case4_bools
# cases_dict['case5'] = case5_bools

# line_colours = ['#000000','#E69F00','#56B4E9', '#FF0E0E', '#30FF0E'] # black, yellow, blue, red, green  
# casenames = ['Eps. MOEA','Int Pen - NodalProp','AOS - NodalProp','Bias Init - NodalProp','ACH - NodalProp']

# alpha_values = [0.5,0.5,0.5,0.5,0.5] # change based on number of cases/visibility

# #nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, nfe_array_1 = hypervolume_computation_all_cases(model_used, cases_dict, artery_problem, sidenum, num_runs, line_colours, alpha_values, casenames, threshold_hv)
# num_fullysat_cases, nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, Uvals_test, Uvals_test_truesat, nfe_array_1 = hypervolume_computation_all_cases(model_used, cases_dict, artery_problem, sidenum, num_runs, line_colours, alpha_values, casenames, threshold_hv)

# ## Mann Whitney U values
# print("For optimization objectives")
# print(Uvals_test)
# print("\n")
# print("For true objectives of fully satisfying designs")
# print(Uvals_test_truesat)

# casenames = ['Eps. MOEA','Int Pen - NodalProp','AOS - NodalProp','Bias Init - NodalProp','ACH - NodalProp']
# plotting_all_cases(num_fullysat_cases, nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, nfe_array_1, line_colours, alpha_values, casenames, threshold_hv)
# #print(nfe_cdf_array)
# #print(hv_dict_truesat_med_cases)



# #### Comparing Simple E-MOEA with Int Pen - Orient, AOS - Orient, Bias Init - Orient and ACH - Orient
# model_used = 2 # 1 = Fibre stiffness, 2 = Truss stiffness, 3 = APDL Beam
# sidenum = 3 # 3x3 node grid
# cases_dict = {}
# artery_problem = False
# num_runs = 30 # number of runs for each case
# threshold_hv = 0.58

# # bools = [int_pen_partcoll, AOS_partcoll, bias_init_partcoll, ACH_partcoll, int_pen_nodalprop, AOS_nodalprop, bias_init_nodalprop, ACH_nodalprop, int_pen_orient, AOS_orient, bias_init_orient, ACH_orient, int_pen_inters, AOS_inters, bias_init_inters, ACH_inters]
# case1_bools = [False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False] # Simple E-MOEA
# case2_bools = [False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False] # Simple E-MOEA
# case3_bools = [False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False] #  
# case4_bools = [False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False] # 
# case5_bools = [False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False] # 

# cases_dict['case1'] = case1_bools
# cases_dict['case2'] = case2_bools
# cases_dict['case3'] = case3_bools
# cases_dict['case4'] = case4_bools
# cases_dict['case5'] = case5_bools

# line_colours = ['#000000','#E69F00','#56B4E9', '#FF0E0E', '#30FF0E'] # black, yellow, blue, red, green  
# casenames = ['Eps. MOEA','Int Pen - Orient','AOS - Orient','Bias Init - Orient','ACH - Orient']

# alpha_values = [0.5,0.5,0.5,0.5,0.5] # change based on number of cases/visibility

# #nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, nfe_array_1 = hypervolume_computation_all_cases(model_used, cases_dict, artery_problem, sidenum, num_runs, line_colours, alpha_values, casenames, threshold_hv)
# num_fullysat_cases, nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, Uvals_test, Uvals_test_truesat, nfe_array_1 = hypervolume_computation_all_cases(model_used, cases_dict, artery_problem, sidenum, num_runs, line_colours, alpha_values, casenames, threshold_hv)

# ## Mann Whitney U values
# print("For optimization objectives")
# print(Uvals_test)
# print("\n")
# print("For true objectives of fully satisfying designs")
# print(Uvals_test_truesat)

# casenames = ['Eps. MOEA','Int Pen - Orient','AOS - Orient','Bias Init - Orient','ACH - Orient']
# plotting_all_cases(num_fullysat_cases, nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, nfe_array_1, line_colours, alpha_values, casenames, threshold_hv)
# #print(nfe_cdf_array)
# #print(hv_dict_truesat_med_cases)



# #### Comparing Simple E-MOEA with Int Pen - Inters, AOS - Inters, Bias Init - Inters and ACH - Inters
# model_used = 2 # 1 = Fibre stiffness, 2 = Truss stiffness, 3 = APDL Beam
# sidenum = 3 # 3x3 node grid
# cases_dict = {}
# artery_problem = False
# num_runs = 30 # number of runs for each case
# threshold_hv = 0.58

# # bools = [int_pen_partcoll, AOS_partcoll, bias_init_partcoll, ACH_partcoll, int_pen_nodalprop, AOS_nodalprop, bias_init_nodalprop, ACH_nodalprop, int_pen_orient, AOS_orient, bias_init_orient, ACH_orient, int_pen_inters, AOS_inters, bias_init_inters, ACH_inters]
# case1_bools = [False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False] # Simple E-MOEA
# case2_bools = [False, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False] # Simple E-MOEA
# case3_bools = [False, False, False, False, False, False, False, False, False, False, False, False, False, True, False, False] #  
# case4_bools = [False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False] # 
# case5_bools = [False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True] # 

# cases_dict['case1'] = case1_bools
# cases_dict['case2'] = case2_bools
# cases_dict['case3'] = case3_bools
# cases_dict['case4'] = case4_bools
# cases_dict['case5'] = case5_bools

# line_colours = ['#000000','#E69F00','#56B4E9', '#FF0E0E', '#30FF0E'] # black, yellow, blue, red, green  
# casenames = ['Eps. MOEA','Int Pen - Inters','AOS - Inters','Bias Init - Inters','ACH - Inters']

# alpha_values = [0.5,0.5,0.5,0.5,0.5] # change based on number of cases/visibility

# #nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, nfe_array_1 = hypervolume_computation_all_cases(model_used, cases_dict, artery_problem, sidenum, num_runs, line_colours, alpha_values, casenames, threshold_hv)
# num_fullysat_cases, nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, Uvals_test, Uvals_test_truesat, nfe_array_1 = hypervolume_computation_all_cases(model_used, cases_dict, artery_problem, sidenum, num_runs, line_colours, alpha_values, casenames, threshold_hv)

# ## Mann Whitney U values
# print("For optimization objectives")
# print(Uvals_test)
# print("\n")
# print("For true objectives of fully satisfying designs")
# print(Uvals_test_truesat)

# casenames = ['Eps. MOEA','Int Pen - Inters','AOS - Inters','Bias Init - Inters','ACH - Inters']
# plotting_all_cases(num_fullysat_cases, nfe_cdf_array, hv_dict_med_cases, hv_dict_1q_cases, hv_dict_3q_cases, hv_dict_true_med_cases, hv_dict_true_1q_cases, hv_dict_true_3q_cases, hv_dict_truesat_med_cases, hv_dict_truesat_1q_cases, hv_dict_truesat_3q_cases, nfe_array_1, line_colours, alpha_values, casenames, threshold_hv)
# #print(nfe_cdf_array)
# #print(hv_dict_truesat_med_cases)




#### Epsilon MOEA runs
# model_used = 2 # 1 = Fibre stiffness, 2 = Truss stiffness, 3 = APDL Beam
# sidenum = 3 # 3x3 node grid
# num_runs = 30 # number of runs for each case
# artery_problem = False

# case_bools = [False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False]
# hypervolume_computation_single_case(model_used, case_bools, artery_problem, sidenum, num_runs, 'Case1_emoea')

stop_time = timeit.default_timer()

print('Run time = ', stop_time - start_time)