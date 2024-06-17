# -*- coding: utf-8 -*-
"""
Hypervolume computation for both satellite problems

@author: roshan94
"""
from pygmo import hypervolume
import csv
import statistics
import numpy as np
import operator as op
from scipy.stats import mannwhitneyu
from functools import reduce
from itertools import combinations
import matplotlib.pyplot as plt
from IPython.core.debugger import set_trace

#### Useful functions and parameter defintions 
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

def get_array_element(array, index):
    return array[index]

#### Determine csv filepath from given case type for one of the satellite problems
def get_csv_filepath_satellite(instrdc_constrained, instrorb_constrained, interinstr_constrained, packeff_constrained, spmass_constrained, instrsyn_constrained, assigning, run_number):
    # instrdc_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # instrorb_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # interinstr_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # packeff_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # spmass_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # instrsyn_constrained = [int_pen, AOS, bias_init, ACH] boolean array
    # assigning = True if assigning problem data is to be read, False if partitioning problem data is to be read
    
    filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\'
    methods = ['Int Pen','AOS','Bias Init','ACH']
    heurs_list = ['Instrdc','Instrorb','Interinstr','Packeff','Spmass','Instrsyn']
    heur_abbrvs_list = ['d','o','i','p','m','s']
    heur_bools = np.vstack((instrdc_constrained, instrorb_constrained, interinstr_constrained, packeff_constrained, spmass_constrained, instrsyn_constrained))
    aos_bools = [x[1] for x in heur_bools]
    
    if (any(aos_bools)):
        filename = 'AOSMOEA_emoea_'
    else:
        filename = 'EpsilonMOEA_emoea_'
        
    if assigning:
        filepath_prob = 'Assigning Problem\\'
        filename_prob = '_assign'
    else:
        filepath_prob = 'Partitioning Problem\\'
        filename_prob = '_partition'
        
    filepath2 = ''
    filename2 = ''
    constr_count = 0
    for i in range(len(heur_bools)):
        constraints = methods[i] + ' - '
        constraints_abbrv = ''
        heur_count = 0
        for j in range(len(heur_bools[0])):
            if heur_bools[j][i]:
                constraints = constraints + heurs_list[j] + '\\'
                constraints_abbrv = constraints_abbrv + heur_abbrvs_list[j]
            else:
                heur_count += 1
            
        if heur_count < len(heur_bools[0]):
            filepath2 = filepath2 + constraints
            filename2 = filename2 + constraints_abbrv + 'con' + str(i) + '_'
        else:
            constr_count += 1
            
    filepath_moea = ''
    if (constr_count == len(heur_bools)):
        filepath_moea = 'Epsilon MOEA\\'
        
    return filepath + filepath_prob + filepath2 + filepath_moea + filename + str(run_number) + filename2 + filename_prob

#### Define NFE array for hypervolume computation (based on number of evaluations in optimization runs)
np.set_printoptions(threshold=np.inf)

# Create array of NFE values at which to compute hypervolume (assumes max function evaluations is 3000)
n_iter_total = 50 # Total number of points in NFE array (1 more than input value to incorporate 0)
n_iter_init = 40 # Number of initial points in NFE array separated by 50 (the rest after that are separated by 100)
nfe_array = np.zeros(n_iter_total+1)
for i in range(n_iter_init):
    nfe_array[i] = 50*i
    
for i in range(n_iter_total - n_iter_init + 1):
    nfe_array[n_iter_init+i] = 50*n_iter_init + 100*i
    
#### Extract Pareto Front and normalization constants data from csv file
def extract_data_from_csv(csv_filepath, assigning, intpen_constr_heur):
    # intpen_constr_heur = [intpen_constr_instrdc, intpen_constr_instrorb, intpen_constr_interinstr, intpen_constr_packeff, intpen_constr_spmass, intpen_constr_instrsyn] boolean array
    with open(csv_filepath,newline='') as csvfile:
        data = [row for row in csv.reader(csvfile)]
                
        num_func_evals_dat = np.zeros(len(data)-1)
        science_pen_dat = np.zeros(len(data)-1)
        cost_pen_dat = np.zeros(len(data)-1)
        science_dat = np.zeros(len(data)-1)
        cost_dat = np.zeros(len(data)-1)
        
        instrdc_scores_dat = np.zeros(len(data)-1)
        instrorb_scores_dat = np.zeros(len(data)-1)
        interinstr_scores_dat = np.zeros(len(data)-1)
        packeff_scores_dat = np.zeros(len(data)-1)
        spmass_scores_dat = np.zeros(len(data)-1)
        instrsyn_scores_dat = np.zeros(len(data)-1)
        
        valid_count = 0
        for x in range(len(data)-1):
            data_float = list(map(float,data[x+1][1:]))
            if (any(np.isnan(np.array(data_float))) or any(np.isinf(np.array(data_float)))):
                continue
            
            num_func_evals_dat[valid_count] = int(data[x+1][1])
            science_pen_dat[valid_count] = int(data[x+1][2])
            cost_pen_dat[valid_count] = int(data[x+1][3])
            science_dat[valid_count] = float(data[x+1][4])
            cost_dat[valid_count] = float(data[x+1][5])
            
            instrdc_scores_dat[valid_count] = float(data[x+1][6])
            instrorb_scores_dat[valid_count] = float(data[x+1][7])
            interinstr_scores_dat[valid_count] = float(data[x+1][8])
            packeff_scores_dat[valid_count] = float(data[x+1][9])
            spmass_scores_dat[valid_count] = float(data[x+1][10])
            instrsyn_scores_dat[valid_count] = float(data[x+1][11])
    
            valid_count += 1
            
    #archs = archs_dat[:valid_count]
    num_func_evals = num_func_evals_dat[:valid_count]
    science_pen = science_pen_dat[:valid_count]
    cost_pen = cost_pen_dat[:valid_count]
    science = science_dat[:valid_count]
    cost = cost_dat[:valid_count]
    
    instrdc_scores = instrdc_scores_dat[:valid_count]
    instrorb_scores = instrorb_scores_dat[:valid_count]
    interinstr_scores = interinstr_scores_dat[:valid_count]
    packeff_scores = packeff_scores_dat[:valid_count]
    spmass_scores = spmass_scores_dat[:valid_count]
    instrsyn_scores = instrsyn_scores_dat[:valid_count]
            
    #print('science')
    #print(science)
    #print('\n')
    #print('cost')
    #print(cost)
    #print('\n')
    
    ## Sort num_fun_evals (and objectives and heuristic scores) in ascending order
    n_func_evals = num_func_evals
    sort_indices = np.argsort(n_func_evals)
    science_pen_sorted = list(science_pen[sort_indices])
    cost_pen_sorted = list(cost_pen[sort_indices])
    cost_sorted = list(cost[sort_indices])
    science_sorted = list(science[sort_indices])
    cost_sorted = list(cost[sort_indices])
    
    instrdc_scores_sorted = list(instrdc_scores[sort_indices])
    instrorb_scores_sorted = list(instrorb_scores[sort_indices])
    interinstr_scores_sorted = list(interinstr_scores[sort_indices])
    packeff_scores_sorted = list(packeff_scores[sort_indices])
    spmass_scores_sorted = list(spmass_scores[sort_indices])
    instrsyn_scores_sorted = list(instrsyn_scores[sort_indices])
    
    #archs_sorted = []
    #for i in range(len(sort_indices)):
        #archs_sorted.append(archs[sort_indices[i]])
    
    heur_scores_sorted = np.vstack((instrdc_scores_sorted, instrorb_scores_sorted, interinstr_scores_sorted, packeff_scores_sorted, spmass_scores_sorted, instrsyn_scores_sorted))
    
    ## Compute only constraint penalized objectives (used for the interior penalty cases)
    heur_obj_weight = 1 # weightage of heuristic penalties wrt objectives
    heur_weights = [[1, 1, 1, 1, 1, 1], [1000, 1000, 1000, 1000, 1000, 1000]] 
    # Row 1 - heur_weights for science, Row 2 - heur_weights for cost
    # heur_weights in each row - [instrdc, instrorb, interinstr, packeff, spmass, instrsyn]
    
    heur_pen = np.zeros((len(instrdc_scores_sorted),2))
    if (any(intpen_constr_heur)):
        heur_index_array = np.arange(len(intpen_constr_heur))
        heur_index_used = [v for i, v in enumerate(heur_index_array) if intpen_constr_heur[i] == True]
     
        for idx in heur_index_used:
            heur_scores_idx = heur_scores_sorted[idx]
            heur_pen_science_idx = [heur_weights[0][idx]*x for x in heur_scores_idx]
            heur_pen_cost_idx = [heur_weights[1][idx]*x for x in heur_scores_idx]
            heur_pen[0] = np.add(heur_pen[0], heur_pen_science_idx/len(heur_index_used))
            heur_pen[1] = np.add(heur_pen[1], heur_pen_cost_idx/len(heur_index_used))
            
    weighted_heur_science_pen = [k*heur_obj_weight for k in heur_pen[0]]
    weighted_heur_cost_pen = [k*heur_obj_weight for k in heur_pen[1]]
    
    science_constr_sorted = list(np.add(science_pen_sorted,weighted_heur_science_pen))
    cost_constr_sorted = list(np.add(cost_pen_sorted,weighted_heur_cost_pen))
    
    ## Determine normalizing objective scores and compute pareto fronts for penalized and true objectives as well as for true objectives of only feasible designs 
    nfe_list_sorted = list(n_func_evals[sort_indices])
    
    #max_func_evals = nfe_list_sorted[-1]
    max_func_evals = 5000 # some runs for some cases run upto 3001 evaluations, which causes hv array length issues

    pareto_front_dict = {}
    #pareto_front_instrdc_dict = {}
    #pareto_front_instrorb_dict = {}
    #pareto_front_interinstr_dict = {}
    #pareto_front_packeff_dict = {}
    #pareto_front_spmass_dict = {}
    #pareto_front_instrsyn_dict = {}
    #pareto_front_archs_dict = {}
    pareto_front_true_dict = {}

    pf_normalize_max_obj1 = []
    pf_normalize_min_obj1 = []
    pf_normalize_max_obj2 = []
    pf_normalize_min_obj2 = []
    pf_true_normalize_max_obj1 = []
    pf_true_normalize_min_obj1 = []
    pf_true_normalize_max_obj2 = []
    pf_true_normalize_min_obj2 = []
    
    count = 0
    pop_size = int(find_last_index(0, nfe_list_sorted))
    nfe_jump_recorded = False

    for nfe_val in nfe_array:
        #print('iter = ' + str(i))
    
        if (nfe_list_sorted[0] == 0):
            if (nfe_val <= 100): # population size set at 100 in java code, but maybe different here due to NaNs
                nfe_index_current = pop_size+1
            else:
                nfe_index_current = find_closest_index(nfe_val, nfe_list_sorted)
        else:
            if (nfe_val <= nfe_list_sorted[0]):
                nfe_index_current = 1
            else:
                nfe_index_current = find_closest_index(nfe_val, nfe_list_sorted)
        
        nfe_array_current = nfe_list_sorted[:nfe_index_current]
        current_population = []
        for j in range(len(nfe_array_current)):
            current_population.append([science_constr_sorted[j], cost_constr_sorted[j]])
            
        #if ("AOS - Orient\\" in csv_filepath) and ("emoea_16" in csv_filepath):
            #set_trace()
        
        current_pareto_front_all = compute_pareto_front(current_population)
        #current_pareto_front = list(set(current_pareto_front_all))
        current_pareto_front = np.unique(current_pareto_front_all, axis=0)
    
        #current_pareto_instrdc_scores = []
        #current_pareto_instrorb_scores = []
        #current_pareto_interinstr_scores = []
        #current_pareto_packeff_scores = []
        #current_pareto_spmass_scores = []
        #current_pareto_instrsyn_scores = []
        #current_pareto_archs = []
        current_pareto_trueobjs = []
        for pareto_arch in current_pareto_front:
            arch_index = science_constr_sorted.index(pareto_arch[0])
            #arch_instrdc_score = get_array_element(instrdc_scores_sorted, arch_index)
            #arch_instrorb_score = get_array_element(instrorb_scores_sorted, arch_index)
            #arch_interinstr_score = get_array_element(interinstr_scores_sorted, arch_index)
            #arch_packeff_score = get_array_element(packeff_scores_sorted, arch_index)
            #arch_spmass_score = get_array_element(spmass_scores_sorted, arch_index)
            #arch_instrsyn_score = get_array_element(instsyn_scores_sorted, arch_index)
            
            #current_pareto_instrdc_scores.append(arch_instrdc_score)
            #current_pareto_instrorb_scores.append(arch_instrorb_score)
            #current_pareto_interinstr_scores.append(arch_interinstr_score)
            #current_pareto_packeff_scores.append(arch_packeff_score)
            #current_pareto_spmass_scores.append(arch_spmass_score)
            #current_pareto_instrsyn_scores.append(arch_instrsyn_score)
            #current_pareto_archs.append(get_array_element(archs_sorted, arch_index))
            
            science_arch, cost_arch = get_true_objectives(science_sorted, cost_sorted, arch_index)
        
            #set_trace()
            
            current_pareto_trueobjs.append([science_arch, cost_arch])
                   
        pareto_front_dict[nfe_val] = current_pareto_front
        #pareto_front_instrdc_dict[nfe_val] = current_pareto_instrdc_scores
        #pareto_front_instrorb_dict[nfe_val] = current_pareto_instrorb_scores
        #pareto_front_interinstr_dict[nfe_val] = current_pareto_interinstr_scores
        #pareto_front_packeff_dict[nfe_val] = current_pareto_packeff_scores
        #pareto_front_spmass_dict[nfe_val] = current_pareto_spmass_scores
        #pareto_front_instrsyn_dict[nfe_val] = current_pareto_instrsyn_scores
        #pareto_front_archs_dict[nfe_val] = current_pareto_archs
        pareto_front_true_dict[nfe_val] = current_pareto_trueobjs
        
        pf_nfeval_obj1 = [row[0] for row in current_pareto_front]
        pf_nfeval_obj2 = [row[1] for row in current_pareto_front]
        
        pf_true_nfeval_obj1 = [row[0] for row in current_pareto_trueobjs]
        pf_true_nfeval_obj2 = [row[1] for row in current_pareto_trueobjs]
        
        pf_normalize_max_obj1.append(np.max(pf_nfeval_obj1))
        pf_normalize_min_obj1.append(np.min(pf_nfeval_obj1))
        pf_normalize_max_obj2.append(np.max(pf_nfeval_obj2))
        pf_normalize_min_obj2.append(np.min(pf_nfeval_obj2))
        
        pf_true_normalize_max_obj1.append(np.max(pf_true_nfeval_obj1))
        pf_true_normalize_min_obj1.append(np.min(pf_true_nfeval_obj1))
        pf_true_normalize_max_obj2.append(np.max(pf_true_nfeval_obj2))
        pf_true_normalize_min_obj2.append(np.min(pf_true_nfeval_obj2))

    ### Computing obj_normalize_fullrun, obj_normalize_afterjump and obj_normalized_true_fullrun using the entire run 
    #obj_normalize_max_fullrun = [np.max(pen_obj1_constr_sorted), np.max(pen_obj2_constr_sorted)]
    #obj_normalize_min_fullrun = [np.min(pen_obj1_constr_sorted), np.min(pen_obj2_constr_sorted)]
    
    #obj_true_normalize_max_fullrun = [np.max(true_obj1_sorted), np.max(true_obj2_sorted)]
    #obj_true_normalize_min_fullrun = [np.min(true_obj1_sorted), np.min(true_obj2_sorted)]

    #obj_normalize_fullrun = [obj_normalize_min_fullrun, obj_normalize_max_fullrun]
    
    #obj_normalize_true_fullrun = [obj_true_normalize_min_fullrun, obj_true_normalize_max_fullrun]
    
    ### Computing obj_normalize_fullrun using the pareto fronts
    obj_normalize_max_fullrun = [np.max(pf_normalize_max_obj1), np.max(pf_normalize_max_obj2)]
    obj_normalize_min_fullrun = [np.min(pf_normalize_min_obj1), np.min(pf_normalize_min_obj2)]  
    
    obj_true_normalize_max_fullrun = [np.max(pf_true_normalize_max_obj1), np.max(pf_true_normalize_max_obj2)]
    obj_true_normalize_min_fullrun = [np.min(pf_true_normalize_min_obj1), np.min(pf_true_normalize_min_obj2)]
    
    obj_normalize_fullrun = [obj_normalize_min_fullrun, obj_normalize_max_fullrun]
    
    obj_normalize_true_fullrun = [obj_true_normalize_min_fullrun, obj_true_normalize_max_fullrun]
    
    return pareto_front_dict, pareto_front_true_dict, obj_normalize_fullrun, obj_normalize_true_fullrun, max_func_evals

#### Compute overall normalization objectives for single case study/all compared case studies (from the complete runs)
def compute_overall_norm_objs(objs_normalization_full, objs_normalization_true):
    # Each input is a dictionary with key as the case study/run string and value as the corresponding 2D array
    
    obj1_max_full_allcases = np.zeros(len(objs_normalization_full))
    obj1_min_full_allcases = np.zeros(len(objs_normalization_full))
    obj2_max_full_allcases = np.zeros(len(objs_normalization_full))
    obj2_min_full_allcases = np.zeros(len(objs_normalization_full))
    obj1_max_true_allcases = np.zeros(len(objs_normalization_true))
    obj1_min_true_allcases = np.zeros(len(objs_normalization_true))
    obj2_max_true_allcases = np.zeros(len(objs_normalization_true))
    obj2_min_true_allcases = np.zeros(len(objs_normalization_true))
    
    i = 0
    for key in objs_normalization_full:
        current_objs_norm_full = objs_normalization_full[key]
        
        obj1_max_full_allcases[i] = current_objs_norm_full[1][0]
        obj2_max_full_allcases[i] = current_objs_norm_full[1][1]
        obj1_min_full_allcases[i] = current_objs_norm_full[0][0]
        obj2_min_full_allcases[i] = current_objs_norm_full[0][1]
        i += 1
        
    i = 0
    for key3 in objs_normalization_true:
        current_objs_norm_true = objs_normalization_true[key3]
        
        obj1_max_true_allcases[i] = current_objs_norm_true[1][0]
        obj2_max_true_allcases[i] = current_objs_norm_true[1][1]
        obj1_min_true_allcases[i] = current_objs_norm_true[0][0]
        obj2_min_true_allcases[i] = current_objs_norm_true[0][1]
        i += 1
        
    obj1_min_full_overall = np.min(obj1_min_full_allcases)
    obj2_min_full_overall = np.min(obj2_min_full_allcases)
    obj1_max_full_overall = np.max(obj1_max_full_allcases)
    obj2_max_full_overall = np.max(obj2_max_full_allcases)
    
    obj1_min_true_overall = np.min(obj1_min_true_allcases)
    obj2_min_true_overall = np.min(obj2_min_true_allcases)
    obj1_max_true_overall = np.max(obj1_max_true_allcases)
    obj2_max_true_overall = np.max(obj2_max_true_allcases)
            
    obj_norm_full_overall = [[obj1_min_full_overall, obj2_min_full_overall], [obj1_max_full_overall, obj2_max_full_overall]]
    obj_norm_true_overall = [[obj1_min_true_overall, obj2_min_true_overall], [obj1_max_true_overall, obj2_max_true_overall]]    
     
    return obj_norm_full_overall, obj_norm_true_overall

#### Compute hypervolume arrays from copmuted pareto fronts and normalization constants
def compute_hv_arrays_from_csv_data(pf_dict, pf_true_dict, obj_norm_full, obj_norm_true_full, max_fun_evals):
    obj_norm_min_full = obj_norm_full[0]
    obj_norm_max_full = obj_norm_full[1]
    obj_norm_true_min_full = obj_norm_true_full[0]
    obj_norm_true_max_full = obj_norm_true_full[1]

    ## Normalize the pareto front objectives and compute the hypervolume
    hypervol_full_dict = []
    hypervol_true_full_dict = []

    for nfe_val in nfe_array:
        #print('iter = ' + str(nfe_val))
    
        current_pareto_front = pf_dict[nfe_val]
        current_true_pareto_front = pf_true_dict[nfe_val]
        current_pf_normalized = []
        current_pf_true_normalized = []
        for pareto_design in current_pareto_front:
            obj1_normalized = (pareto_design[0] - obj_norm_min_full[0])/(obj_norm_max_full[0] - obj_norm_min_full[0])
            obj2_normalized = (pareto_design[1] - obj_norm_min_full[1])/(obj_norm_max_full[1] - obj_norm_min_full[1])
            current_pf_normalized.append([obj1_normalized, obj2_normalized])
                    
        for pareto_design_true in current_true_pareto_front:
            obj1_true_normalized = (obj_norm_true_max_full[0] - pareto_design_true[0])/(obj_norm_true_max_full[0] - obj_norm_true_min_full[0])
            obj2_true_normalized = (pareto_design_true[1] - obj_norm_true_min_full[1])/(obj_norm_true_max_full[1] - obj_norm_true_min_full[1])
            current_pf_true_normalized.append([obj1_true_normalized, obj2_true_normalized])
            
        current_hv = compute_hv(current_pf_normalized)
        hypervol_full_dict.append([nfe_val, current_hv])
        
        current_hv_true = compute_hv(current_pf_true_normalized)
        hypervol_true_full_dict.append([nfe_val, current_hv_true])
        
    return hypervol_full_dict, hypervol_true_full_dict

### Compute array of NFE values for reaching threshold hypervolume for different runs of a particular case
def compute_nfe_hypervolume_attained(hv_dict):
    hv_threshold = 0.75 # Threshold HV value to reach, user parameter
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
def plot_fraction_hypervolume_attained(nfe_hv_attained_dict, nfe_array, colour_array, casename_array, savefig_name):
    fig1 = plt.figure(1)
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
    plt.ylabel('Fraction of runs HV $\geq$ 0.75',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.20), ncol=3, borderaxespad=0 ,prop={"size":12})
    plt.show()
    fig1.savefig('frac_hv_attained_' + savefig_name + '.png', format='png')
    
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

def plot_hypervolume_stats(hv_median_case, hv_1q_case, hv_3q_case, nfe_array, savefig_name):
    fig1 = plt.figure(1)
    plt.plot(nfe_array,hv_median_case, 'b-', label='Median')
    plt.plot(nfe_array,hv_1q_case, 'r-', label='1st Quartile')
    plt.plot(nfe_array,hv_3q_case, 'g-', label='3rd Quartile')
    plt.xlabel('Number of Function Evaluations')
    plt.ylabel('Hypervolume')
    plt.title('Averaged Hypervolume vs NFE')
    plt.legend(loc='lower right')
    plt.show()
    #fig1.savefig('HV_plot_averaged_' + savefig_name + '.png')
    
def plot_hypervolume_stats_allcases(hv_median_dict, hv_1q_dict, hv_3q_dict, nfe_array, colour_array, alpha_array, casename_array, plot_title, savefig_name):
    fig1 = plt.figure(1)
    number_cases = len(hv_median_dict)
    #print('n_cases')
    #print(number_cases)
    for i in range(number_cases):
        #print(print(marker_array[i]+'*'))
        plt.plot(nfe_array, hv_median_dict['case'+str(i)], '-', color=colour_array[i], label=casename_array[i])
        #plt.fill_between(nfe_array, hv_1q_dict['case'+str(i)], hv_3q_dict['case'+str(i)], color=colour_array[i], alpha=alpha_array[i])
        
        plt.plot(nfe_array, hv_1q_dict['case'+str(i)], '--', color=colour_array[i])#, label=casename_array[i]+' 1st Quartile')
        plt.plot(nfe_array, hv_3q_dict['case'+str(i)], '--', color=colour_array[i])#, label=casename_array[i]+' 3rd Quartile')
    plt.xlabel('Number of Function Evaluations',fontsize=14)
    plt.ylabel('Hypervolume',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.title(plot_title)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.20), ncol=3, borderaxespad=0, prop={"size":12})
    plt.show()
    fig1.savefig('HV_plot_averaged_' + savefig_name + '.png', format='png')
    
def compute_mann_whitney_Uvals(hv_med_dict_allcases, nfe_array): # Wilcoxon Rank Sum Test
    # hv_med_array_allcases is a dictionary of length = number of cases
    n_samples = 20
    linspace_samples_array = np.linspace(0,1,n_samples)
    nfe_samples_array = np.floor(np.multiply(linspace_samples_array, nfe_array[-1]))
    nfe_samples_indices_array = np.zeros(len(nfe_samples_array))
    for i in range(len(nfe_samples_array)):
        nfe_samples_indices_array[i] = find_closest_index(nfe_samples_array[i], nfe_array)
        
    hv_med_samples_allcases = {}
    for j in range(len(hv_med_dict_allcases)):
        hv_med_case = hv_med_dict_allcases['case'+str(j)]
        hv_med_samples_case = np.zeros(len(hv_med_case))
        for k in range(len(hv_med_case)):
            hv_med_samples_case[k] = hv_med_case[nfe_samples_indices_array[k]]
        hv_med_samples_allcases['case'+str(j)] = hv_med_samples_case
        
    cases_array = np.arange(len(hv_med_dict_allcases))
    case_combinations = list(combinations(cases_array,2))
    
    U_test_cases = {}
    
    for k in range(len(case_combinations)):
        case_string_x = 'case' + str(case_combinations[k][0])
        case_string_y = 'case' + str(case_combinations[k][1])
        
        U1, p = mannwhitneyu(hv_med_samples_allcases[case_string_x], hv_med_samples_allcases[case_string_y])
        U2 = len(hv_med_samples_allcases[case_string_x])*len(hv_med_samples_allcases[case_string_y]) - U1
        
        U_test = np.min(np.array([U1, U2]))
        dict_key = case_string_x + ' and ' + case_string_y
        
        U_test_cases[dict_key] = [U_test, p]
        
    return U_test_cases

#### Define functions to compute and plot hypervolume for single case and all cases
def hypervolume_computation_single_case(case_booleans, prob_assigning, run_nums, case_name):
    ## Computing the pareto fronts and normalization objectives for each run
    obj_norm_allruns = {}
    obj_norm_true_allruns = {}
    pf_allruns = {}
    pf_true_allruns = {}
    max_f_evals_allruns = np.zeros(run_nums)
    for i in range(run_nums):
        print('Computing Pareto Fronts for run ' + str(i))
        current_csvpath = get_csv_filepath_satellite(case_booleans[:4], case_booleans[4:8], case_booleans[8:12], case_booleans[12:16], case_booleans[16:20], case_booleans[20:24], prob_assigning, i)
        heur_intpen_constr = [case_booleans[0], case_booleans[4], case_booleans[8], case_booleans[12]]
        pf_dict_i, pf_true_dict_i, obj_norm_full_i, obj_norm_true_i, max_fun_evals_i = extract_data_from_csv(current_csvpath, prob_assigning, heur_intpen_constr)
        pf_allruns['run'+str(i)] = pf_dict_i
        pf_true_allruns['run'+str(i)] = pf_true_dict_i
        obj_norm_allruns['run'+str(i)] = obj_norm_full_i
        obj_norm_true_allruns['run'+str(i)] = obj_norm_true_i
        max_f_evals_allruns[i] = max_fun_evals_i
    
    #print('pf_allruns')
    #print(pf_allruns)
    #print('\n')
    #print('pf_true_allruns')
    #print(pf_true_allruns)
    #print('\n')
    ## Use computed normalization objectives and find the overall normalization objectives across all runs
    print('Computing overall normalization constants')
    norm_objs_full_overall, norm_objs_true_overall = compute_overall_norm_objs(obj_norm_allruns, obj_norm_true_allruns)
    #print('norm_objs_full_overall')
    #print(norm_objs_full_overall)
    #print('\n')
    #print('norm_objs_true_overall')
    #print(norm_objs_true_overall)
    #print('\n')
    
    ## Compute Hypervolume values for each run
    hv_dict_allruns = {}
    hv_dict_true_allruns = {}
    for j in range(run_nums):
        print('Computing hypervolume values for run ' + str(j))
        hv_dict_j, hv_dict_true_j = compute_hv_arrays_from_csv_data(pf_allruns['run'+str(j)], pf_true_allruns['run'+str(j)], norm_objs_full_overall, norm_objs_true_overall, max_f_evals_allruns[j])
        hv_dict_allruns['run'+str(j)] = hv_dict_j
        hv_dict_true_allruns['run'+str(j)] = hv_dict_true_j
        
    #print('hv_dict_allruns')
    #print(hv_dict_allruns)
    #print('hv_dict_true_allruns')
    #print(hv_dict_true_allruns)
        
    ## Plotting
    print('Plotting')
    hv_median_all, hv_1q_all, hv_3q_all, nfe_array = compute_hypervolume_stats(hv_dict_allruns)
    plot_hypervolume_stats(hv_median_all, hv_1q_all, hv_3q_all, nfe_array, case_name+'_full')

    ## Plot HVs for hv_afterjump

    ## Plot HVs for true objectives
    hv_true_median_all, hv_true_1q_all, hv_true_3q_all, nfe_array_true = compute_hypervolume_stats(hv_dict_true_allruns)
    plot_hypervolume_stats(hv_true_median_all, hv_true_1q_all, hv_true_3q_all, nfe_array_true, case_name+'_true')
    
    
def hypervolume_computation_all_cases(choice_model, case_bools_dict, prob_assigning, run_nums, marker_colours, alpha_vals, case_names):
    num_cases = len(case_bools_dict) # number of cases to compare 

    ## Computing the pareto fronts and normalization objectives for each run in each case
    pf_allcases = {}
    pf_true_allcases = {}
    obj_norm_allcasesandruns = {}
    obj_norm_true_allcasesandruns = {}
    max_f_evals_allcases = {}
    for i in range(num_cases):
        print('Computing Pareto Fronts for runs in Case '+str(i))
        current_case_bools = case_bools_dict['case'+str(i+1)]
        #set_trace()
        pf_allruns_i = {}
        pf_true_allruns_i = {}
        max_f_evals_allruns = np.zeros(run_nums)
        for j in range(run_nums):
            print('Run '+str(j))
            current_csvpath = get_csv_filepath_satellite(current_case_bools[:4], current_case_bools[4:8], current_case_bools[8:12], current_case_bools[12:16], current_case_bools[16:20], current_case_bools[20:24], prob_assigning, j)
            #set_trace()
            heur_intpen_constr = [current_case_bools[0], current_case_bools[4], current_case_bools[8], current_case_bools[12]]
            pf_dict_j, pf_true_dict_j, obj_norm_full_j, obj_norm_true_j, max_fun_evals_j = extract_data_from_csv(current_csvpath, prob_assigning, heur_intpen_constr)
            pf_allruns_i['run'+str(j)] = pf_dict_j
            pf_true_allruns_i['run'+str(j)] = pf_true_dict_j
            obj_norm_allcasesandruns['case'+str(i+1)+'run'+str(j)] = obj_norm_full_j
            obj_norm_true_allcasesandruns['case'+str(i+1)+'run'+str(j)] = obj_norm_true_j
            max_f_evals_allruns[j] = max_fun_evals_j
        pf_allcases['case'+str(i+1)] = pf_allruns_i
        pf_true_allcases['case'+str(i+1)] = pf_true_allruns_i
        max_f_evals_allcases['case'+str(i+1)] = max_f_evals_allruns
    
    ## Use computed normalization objectives and find the overall normalization objectives across all runs and cases
    print('Computing overall normalization constants')
    norm_objs_full_overall, norm_objs_true_overall = compute_overall_norm_objs(obj_norm_allcasesandruns, obj_norm_true_allcasesandruns)
    
    #set_trace()
    
    ## Compute Hypervolume values for each run in each case
    hv_dict_median_allcases = {}
    hv_dict_1q_allcases = {}
    hv_dict_3q_allcases = {}
    hv_dict_true_median_allcases = {}
    hv_dict_true_1q_allcases = {}
    hv_dict_true_3q_allcases = {}
    nfe_array_hv_attained_dict = {}
    for i in range(num_cases):
        print('Computing hypervolume values for runs in Case '+str(i))
        pfs_case_i = pf_allcases['case'+str(i+1)]
        pfs_true_case_i = pf_true_allcases['case'+str(i+1)]
        max_func_evals_i = max_f_evals_allcases['case'+str(i+1)]
        hv_dict_allruns = {}
        hv_dict_true_allruns = {}
        for j in range(run_nums):
            print('Run '+str(j))
            hv_dict_j, hv_dict_true_j = compute_hv_arrays_from_csv_data(pfs_case_i['run'+str(j)], pfs_true_case_i['run'+str(j)], norm_objs_full_overall, norm_objs_true_overall, max_func_evals_i[j])
            hv_dict_allruns['run'+str(j)] = hv_dict_j
            hv_dict_true_allruns['run'+str(j)] = hv_dict_true_j
            
        print('Computing array of NFE for attaining threshold hypervolume')
        nfe_hv_attained_case = compute_nfe_hypervolume_attained(hv_dict_allruns)
        nfe_array_hv_attained_dict['case'+str(i)] = nfe_hv_attained_case
                
        print('Computing hypervolume stats')
        hv_med_i, hv_1q_i, hv_3q_i, nfe_array_i = compute_hypervolume_stats(hv_dict_allruns)
        #hv_med_aj_i, hv_1q_aj_i, hv_3q_aj_i, nfe_array_aj_i = compute_hypervolume_stats(hv_dict_aj_allruns)
        hv_med_true_i, hv_1q_true_i, hv_3q_true_i, nfe_array_true_i = compute_hypervolume_stats(hv_dict_true_allruns)
                
        hv_dict_median_allcases['case'+str(i)] = hv_med_i
        hv_dict_1q_allcases['case'+str(i)] = hv_1q_i
        hv_dict_3q_allcases['case'+str(i)] = hv_3q_i
        hv_dict_true_median_allcases['case'+str(i)] = hv_med_true_i
        hv_dict_true_1q_allcases['case'+str(i)] = hv_1q_true_i
        hv_dict_true_3q_allcases['case'+str(i)] = hv_3q_true_i
          
    return nfe_array_hv_attained_dict, hv_dict_median_allcases, hv_dict_1q_allcases, hv_dict_3q_allcases, hv_dict_true_median_allcases, hv_dict_true_1q_allcases, hv_dict_true_3q_allcases, nfe_array_i
    
def plotting_all_cases(nfe_hv_attained_dict, hv_dict_med_allcases, hv_dict_1stq_allcases, hv_dict_3rdq_allcases, hv_dict_true_med_allcases, hv_dict_true_1stq_allcases, hv_dict_true_3rdq_allcases, nfe_array0, mark_colors, alphas, names_cases):
    print('Plotting')
    plot_fraction_hypervolume_attained(nfe_hv_attained_dict, nfe_array0, mark_colors, names_cases, 'allcases_full')
    plot_hypervolume_stats_allcases(hv_dict_med_allcases, hv_dict_1stq_allcases, hv_dict_3rdq_allcases, nfe_array, mark_colors, alphas, names_cases, 'Hypervolume of Penalized Objectives', 'allcases_full')
    plot_hypervolume_stats_allcases(hv_dict_true_med_allcases, hv_dict_true_1stq_allcases, hv_dict_true_3rdq_allcases, nfe_array, mark_colors, alphas, names_cases, 'Hypervolume of True Objectives', 'allcases_true')
    
