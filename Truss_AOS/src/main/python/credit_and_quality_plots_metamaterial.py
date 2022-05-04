# -*- coding: utf-8 -*-
"""
Plot credits and quality values for each heuristic averaged across all runs of 
the metamaterial problems

@author: roshan94
"""
import numpy as np
import csv
import statistics
import matplotlib.pyplot as plt

artery_problem = True
model = 1 # 0 -> Fibre model, 1 -> Truss model, 2 -> Beam model

### Read the credit and quality csv files
file_loc = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\'

aos_heur_bools = [True, True, True, True] # [part_coll, nodal_prop, orient, inters]

credit_assign = 2 # 0 -> offspring parent dominance, 1 -> set improvement dominance, 2 -> set contribution dominance

#run_num = 1 # run number of results to read
num_runs = 30

heurs_list = ['PartColl','NodalProp','Orient','Inters']
heur_abbrvs_list = ['p','n','o','i']

filename = 'AOSMOEA_emoea_'

filename_prob = '_prob2'
filepath_prob = 'Truss Problem\\'
if artery_problem:
    filename_prob = '_artery'
    filepath_prob = 'Artery Problem\\'
    
filename_model = '_fibre'
filepath_model = 'Fibre Model\\'
if model == 1:
    filename_model = '_truss'
    filepath_model = 'Truss Model\\'
elif model == 2:
    filename_model = '_beam'
    filepath_model = 'Beam Model\\'
    
heurs = ''
heurs_path = 'AOS - '
for i in range(len(aos_heur_bools)):
    if aos_heur_bools[i]:
        heurs = heurs + heur_abbrvs_list[i]
        heurs_path = heurs_path + heurs_list[i] 

if all(aos_heur_bools):
    operator_strings = ['AddDiagonalMember+BitFlip','AddMember+BitFlip','ImproveOrientation2+BitFlip','RemoveIntersection2+BitFlip','OnePointCrossover+BitFlip']
else:
    operator_strings = ['ImproveOrientation2+BitFlip','RemoveIntersection2+BitFlip','OnePointCrossover+BitFlip']

filepath_cred = 'offspring parent dominance\\'
if credit_assign == 1:
    filepath_cred = 'set improvement dominance\\'
elif credit_assign == 2:
    filepath_cred = 'set contribution dominance\\'

pop_size = 100
max_func_eval = 6000

nfe_array = np.linspace(pop_size+2, max_func_eval, (int((max_func_eval-(pop_size+2))/2)+1))
nfe_array = nfe_array.astype(int)

def get_csv_rows(filepath):
    rows_dict = {}
    
    with open(filepath, newline='') as csvfile:
        file_reader = csv.reader(csvfile)
        
        ind_row = 0
        for row in file_reader:
            rows_dict[ind_row] = row
            ind_row = ind_row + 1
            
    return rows_dict

def get_value_dicts_run(filename_val, op_strings, nfe_arr):
    # value can be credits or qualities, assembled for each operator
    val_rows_filename = get_csv_rows(filename_val)
    vals_run = {}
    for j in range(len(op_strings)):
        val_array_operator = np.zeros((len(nfe_arr)))
        for k in range(len(val_rows_filename)):
            if (val_rows_filename[k][0] == op_strings[j]):
                val_operator = val_rows_filename[k]
                nfe_operator = val_rows_filename[k-1]
                for m in range(len(nfe_arr)):
                    if (str(nfe_arr[m]) in nfe_operator):
                        idx = nfe_operator.index(str(nfe_arr[m]))
                        val_array_operator[m] = float(val_operator[idx])
                    else:
                        val_array_operator[m] = 0
            else:
                continue
        vals_run[op_strings[j]] = val_array_operator
    return vals_run

def plot_values(val_mean_dict, val_std_dict, op_strings, labels, nfe_arr, ylab):
    fig = plt.figure()
    for i in range(len(op_strings)):
        val_array_mean = val_mean_dict[op_strings[i]]
        val_array_std = val_std_dict[op_strings[i]]
        plt.plot(nfe_arr, val_array_mean, label=labels[i])
        val_plus_err = np.add(val_array_mean, np.multiply(val_array_std, 2))
        val_minus_err = np.subtract(val_array_mean, np.multiply(val_array_std, 2))
        plt.fill_between(nfe_arr, val_minus_err, val_plus_err, alpha=0.5)
    plt.xlabel('NFE')
    plt.ylabel(ylab)
    plt.legend(loc="best")
    plt.show()

cred_vals = {}
qual_vals = {}

### Credit arrays for each operator for all runs
for i in range(num_runs):
    cred_filename = file_loc + filepath_prob + 'Constant Radii\\' + filepath_model + heurs_path + '\\' + filepath_cred + 'emoea_' + str(i) + heurs + 'con1_' + filename_prob + filename_model + '_credit.csv'
    cred_vals_run = get_value_dicts_run(cred_filename, operator_strings, nfe_array)
    cred_vals[i] = cred_vals_run
            
### Quality arrays for each operator for all runs        
for i in range(num_runs):
    qual_filename = file_loc + filepath_prob + 'Constant Radii\\' + filepath_model + heurs_path + '\\' + filepath_cred + 'emoea_' + str(i) + heurs + 'con1_' + filename_prob + filename_model + '_qual.csv'
    qual_vals_run = get_value_dicts_run(qual_filename, operator_strings, nfe_array)
    qual_vals[i] = qual_vals_run

### Average credits and quality values over each run
creds_arrays_mean = {}
creds_arrays_std = {}
quals_arrays_mean = {}
quals_arrays_std = {}

for i in range(len(operator_strings)):
    creds_operator_mean = np.zeros((len(nfe_array)))
    creds_operator_std = np.zeros((len(nfe_array)))
    
    quals_operator_mean = np.zeros((len(nfe_array)))
    quals_operator_std = np.zeros((len(nfe_array)))
    for j in range(len(nfe_array)):
        cred_operator_nfe = np.zeros((num_runs))
        qual_operator_nfe = np.zeros((num_runs))
        for k in range(num_runs):
            creds_run = cred_vals[k]
            quals_run = qual_vals[k]
            
            creds_run_operator = creds_run[operator_strings[i]]
            quals_run_operator = quals_run[operator_strings[i]]
            
            cred_operator_nfe[k] = creds_run_operator[j]
            qual_operator_nfe[k] = quals_run_operator[j]
        
        creds_operator_mean[j] = statistics.mean(cred_operator_nfe)
        creds_operator_std[j] = statistics.stdev(cred_operator_nfe)
        
        quals_operator_mean[j] = statistics.mean(qual_operator_nfe)
        quals_operator_std[j] = statistics.stdev(qual_operator_nfe)
            
    creds_arrays_mean[operator_strings[i]] = creds_operator_mean
    creds_arrays_std[operator_strings[i]] = creds_operator_std
    
    quals_arrays_mean[operator_strings[i]] = quals_operator_mean
    quals_arrays_std[operator_strings[i]] = quals_operator_std
            
### Plot credits and quality values
if all(aos_heur_bools):
    heur_labels = ['PartColl','NodalProp','Orient','Inters','X+M']
else:
    heur_labels = ['Orient','Inters','X+M']

# Credit plots
plot_values(creds_arrays_mean, creds_arrays_std, operator_strings, heur_labels, nfe_array, 'Credit')

# Quality plots
plot_values(quals_arrays_mean, quals_arrays_std, operator_strings, heur_labels, nfe_array, 'Quality')
    
    
        
        
    
                
            
    
    

