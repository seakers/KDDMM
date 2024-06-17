# -*- coding: utf-8 -*-
"""
Biased sampling screening study

@author: roshan94
"""
from pygmo import hypervolume
import csv
import numpy as np
from IPython.core.debugger import set_trace

### Useful functions
def get_csv_filepath_material(bias_init_enforced, artery_problem, model_choice, run_number):
    # bias_init_enforced = [partcoll, nodalprop, orient, inters] boolean array
    # artery_problem = True if artery problem data is to be read, False if truss problem data is to be read
    # model_choice = 1 - fibre stiffness, 2 - truss stiffness, 3 - ANSYS APDL beam model
    
    filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\' # for office system
    #filepath = 'C:\\Users\\rosha\\Documents\\SEAK Lab Github\\KD3M3\\result\\' # for home system
    heurs_list = ['Partial Collapsibility','Nodal Properties','Orientation','Intersection']
    heur_abbrvs_list = ['pc','np','or','is']
    
    if artery_problem:
        filepath_prob = 'Artery Problem\\'
        filename_prob = '_artery'
    else:
        filepath_prob = 'Truss Problem\\'
        filename_prob = '_prob2'
        
    filename = 'random_biased_sampling_index_data'
        
    filepath2 = ''
    filename2 = ''
    constraints = ''
    constraints_abbrv = ''
    constr_count = 0
    for i in range(len(bias_init_enforced)):
        if bias_init_enforced[i]:
            constraints = constraints + heurs_list[i]
            constraints_abbrv = constraints_abbrv + heur_abbrvs_list[i]
        else:
            constr_count += 1
    
    filepath_rand = ''        
    if constr_count < len(bias_init_enforced):
        filepath2 = filepath2 + constraints + '\\'
        filename2 = filename2 + '_' + constraints_abbrv 
    else:
        filepath_rand = 'random\\'
            
    if model_choice == 1:
        filepath_model = 'Fiber Model\\'
        filename_model = '_fibre'
    elif model_choice == 2:
        filepath_model = 'Truss Model\\'
        filename_model = '_truss'
    elif model_choice == 3:
        filepath_model = 'Beam Model\\'
        filename_model = '_beam'
        
    return filepath + 'Biased Sampling Index Data\\' + filepath_prob + filepath_model + filepath2 + filepath_rand + filename + filename_prob + filename_model + filename2 + str(run_number) + '.csv' 

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
    
    #pen_weight = np.zeros((3))
    #pen_weight[0] = np.mean(feas_viol_array) # feas_viol_arrray is already in [0, 1]
    #pen_weight[1] = np.mean(conn_viol_array) # conn_viol_arrray is already in [0, 1]
    #pen_weight[2] = np.mean(stiffrat_array) # stiffrat_array is roughly between [0, 1] for truss problem
    
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

def find_last_index(val,search_list):
    if val in search_list:
        idx = len(search_list) - search_list[::-1].index(val) - 1
    else:
        idx = 0
    return idx

def compute_hv(population):
    array_archs = np.zeros((len(population), 2))
    for i in range(len(population)):
        array_archs[i] = population[i]
    hv_object = hypervolume(array_archs)
    hv = hv_object.compute([1.1,1.1])/1.1**2
    return hv

def extract_data_from_csv(csv_filepath, artery_problem):
    # intpen_constr_heur = [intpen_constr_partcoll, intpen_constr_nodalprop, intpen_constr_orient, intpen_constr_inters] boolean array
    with open(csv_filepath,newline='') as csvfile:
        data = [row for row in csv.reader(csvfile)]
                
        norm_obj1_dat = np.zeros(len(data)-1)
        norm_obj2_dat = np.zeros(len(data)-1)
        
        feas_scores_dat = np.zeros(len(data)-1)
        conn_scores_dat = np.zeros(len(data)-1)
        if (not artery_problem):
            stiffrat_vals_dat = np.zeros(len(data)-1)
        
        valid_count = 0
        for x in range(len(data)-1):
            data_float = list(map(float,data[x+1][1:-1]))
            if (any(np.isnan(np.array(data_float))) or any(np.isinf(np.array(data_float)))):
                continue
            
            norm_obj1_dat[valid_count] = float(data[x+1][1])
            norm_obj2_dat[valid_count] = float(data[x+1][2])
            
            feas_scores_dat[valid_count] = float(data[x+1][3])
            conn_scores_dat[valid_count] = float(data[x+1][4])
            if (not artery_problem):
                if (float(data[x+1][5]) > 50):
                    stiffrat_vals_dat[valid_count] = 10
                else:
                    stiffrat_vals_dat[valid_count] = float(data[x+1][5])
    
            valid_count += 1
            
    #designs = designs_dat[:valid_count]
    norm_obj1 = norm_obj1_dat[:valid_count]
    norm_obj2 = norm_obj2_dat[:valid_count]
    
    feas_scores =  feas_scores_dat[:valid_count]
    conn_scores = conn_scores_dat[:valid_count]
    if (not artery_problem):
        stiffrat_vals = stiffrat_vals_dat[:valid_count]
    else:
        stiffrat_vals = np.zeros((len(feas_scores)))
    
    ## Get population
    norm_objs_init = np.column_stack((norm_obj1, norm_obj2))
    pen_objs_init = get_pen_objs(norm_objs_init, feas_scores, conn_scores, stiffrat_vals)
    
    return pen_objs_init

def find_norm_vals(pop):
    pen_obj1_max = np.amax([x[0] for x in pop])
    pen_obj1_min = np.amin([x[0] for x in pop])
    pen_obj2_max = np.amax([x[1] for x in pop])
    pen_obj2_min = np.amin([x[1] for x in pop])
    
    return pen_obj1_max, pen_obj1_min, pen_obj2_max, pen_obj2_min

#### OPERATION
artery_prob = True # True for artery problem, False for truss problem
model_used = 2 # 1 = Fibre stiffness, 2 = Truss stiffness, 3 = APDL Beam

eps_moea_enforced = [False, False, False, False]
orient_enforced = [False, False, True, False]

bias_sampling_cases = np.vstack((eps_moea_enforced, orient_enforced))

n_runs = 10
n_cases = bias_sampling_cases.shape[0]

## Extract and store initial populations
init_pop_allcases = {}
for i in range(n_cases):
    init_pop_allruns = {}
    for j in range(n_runs):
        filepath_run = get_csv_filepath_material(bias_sampling_cases[i], artery_prob, model_used, j)
        init_pop_run = extract_data_from_csv(filepath_run, artery_prob)
        init_pop_allruns['run'+str(j)] = init_pop_run
    init_pop_allcases['case'+str(i)] = init_pop_allruns
    
## Find overall normalization constants
pen_obj1_max_all = np.zeros((n_cases, n_runs))
pen_obj1_min_all = np.zeros((n_cases, n_runs))
pen_obj2_max_all = np.zeros((n_cases, n_runs))
pen_obj2_min_all = np.zeros((n_cases, n_runs))

for i in range(n_cases):
    init_pop_caseruns = init_pop_allcases['case'+str(i)]
    for j in range(n_runs):
        init_pop_run = init_pop_caseruns['run'+str(j)]
        penobj1_max, penobj1_min, penobj2_max, penobj2_min = find_norm_vals(init_pop_run)
        pen_obj1_max_all[i][j] = penobj1_max
        pen_obj1_min_all[i][j] = penobj1_min
        pen_obj2_max_all[i][j] = penobj2_max
        pen_obj2_min_all[i][j] = penobj2_min
    
pen_obj1_norm_max = np.amax(pen_obj1_max_all)
pen_obj1_norm_min = np.amin(pen_obj1_min_all)
pen_obj2_norm_max = np.amax(pen_obj2_max_all)
pen_obj2_norm_min = np.amin(pen_obj2_min_all)

#pen_obj1_norm_max = np.zeros((n_runs))
#pen_obj1_norm_min = np.zeros((n_runs))
#pen_obj2_norm_max = np.zeros((n_runs))
#pen_obj2_norm_min = np.zeros((n_runs))

#for i in range(n_runs):
    #pen_obj1_norm_max[i] = np.amax(pen_obj1_max_all[:,i])
    #pen_obj1_norm_min[i] = np.amin(pen_obj1_min_all[:,i])
    #pen_obj2_norm_max[i] = np.amax(pen_obj2_max_all[:,i])
    #pen_obj2_norm_min[i] = np.amin(pen_obj2_min_all[:,i])

## Compute hypervolumes
hv_cases = np.zeros((n_runs, n_cases))
for i in range(n_cases):
    init_pop_caseruns = init_pop_allcases['case'+str(i)]
    for j in range(n_runs):
        init_pop_run = init_pop_caseruns['run'+str(j)]
        pen_obj1_run = [x[0] for x in init_pop_run]
        pen_obj2_run = [x[1] for x in init_pop_run]
        
        pen_obj1_norm_run = (pen_obj1_run - pen_obj1_norm_min)/(pen_obj1_norm_max - pen_obj1_norm_min)
        pen_obj2_norm_run = (pen_obj2_run - pen_obj2_norm_min)/(pen_obj2_norm_max - pen_obj2_norm_min)
        
        #pen_obj1_norm_run = (pen_obj1_run - pen_obj1_norm_min[j])/(pen_obj1_norm_max[j] - pen_obj1_norm_min[j])
        #pen_obj2_norm_run = (pen_obj2_run - pen_obj2_norm_min[j])/(pen_obj2_norm_max[j] - pen_obj2_norm_min[j])
        
        pop_norm = np.column_stack((pen_obj1_norm_run, pen_obj2_norm_run))
        hv_cases[j,i] = compute_hv(pop_norm)
        
## Compute indices
indices = np.zeros((n_runs, n_cases-1))
for i in range(n_cases-1):
    indices = hv_cases[:,i+1] - hv_cases[:,0]

print('Orientation Impact Index: ' + str(np.mean(indices)) + ' +/- ' + str(np.std(indices)))

### Computing minumum percentile for positive indices 
percentile_vals = np.linspace(1, 100, 100)
if hv_cases.shape[1] == 2:
    n_perctile_heurs = [0]
else:    
    n_perctile_heurs = np.zeros((indices.shape[1]))
    
for i in range(len(n_perctile_heurs)):
    if hv_cases.shape[1] == 2:
        I_current_heur = indices
    else:
        I_current_heur = indices[:,i]
    for j in range(len(percentile_vals)):
        pctile = np.percentile(I_current_heur, percentile_vals[j], method='interpolated_inverted_cdf')
        if pctile > 0:
            n_perctile_heurs[i] = percentile_vals[j]
            break
        if j == len(percentile_vals)-1:
            n_perctile_heurs[i] = percentile_vals[j]
