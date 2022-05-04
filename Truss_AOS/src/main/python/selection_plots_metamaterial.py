# -*- coding: utf-8 -*-
"""
Plot selection histories for each heuristic for the metamaterial problems

@author: roshan94
"""
import numpy as np
import csv
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

nfe_array = np.linspace(pop_size, max_func_eval, (int((max_func_eval-pop_size)/2)+1))
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

def get_selection_histories_run(filename_sel, op_strings):
    sel_rows_filename = get_csv_rows(filename_sel)
    op_sels = {}
    
    for i in range(len(op_strings)):
        op_sel_nfes = []
        for j in range(len(sel_rows_filename)):
            if (sel_rows_filename[j][1] == op_strings[i]):
                op_sel_nfes.append(sel_rows_filename[j][0])
        op_sels[op_strings[i]] = op_sel_nfes
    
    return op_sels
    
### Get operator selection histories for each run
sel_hists_runs = {}
for i in range(num_runs):
    sel_filename = file_loc + filepath_prob + 'Constant Radii\\' + filepath_model + heurs_path + '\\' + filepath_cred + 'emoea_' + str(i) + heurs + 'con1_' + filename_prob + filename_model + '_hist.csv'
    ops_hist_run = get_selection_histories_run(sel_filename, operator_strings)
    sel_hists_runs[i] = ops_hist_run
        
### Get selection frequencies for each heuristic across all runs
sel_freq = {}
for i in range(len(operator_strings)):
    sel_freq_op = np.zeros((len(nfe_array)))
    for j in range(len(nfe_array)):
        num_sel_op = 0
        for k in range(num_runs):
            sel_hist_run = sel_hists_runs[k]
            sel_hist_op_run = sel_hist_run[operator_strings[i]]
            if (str(nfe_array[j]) in sel_hist_op_run):
                num_sel_op = num_sel_op + 1
        sel_freq_op[j] = num_sel_op/num_runs
    sel_freq[operator_strings[i]] = sel_freq_op
        
### Plot selection frequencies as stacked bar graph
if all(aos_heur_bools):
    heur_labels = ['PartColl','NodalProp','Orient','Inters','X+M']
else:
    heur_labels = ['Orient','Inters','X+M']
    
fig = plt.figure()
sel_freq_sum = sel_freq[operator_strings[0]]
plt.bar(nfe_array, sel_freq_sum)
for i in range(1, len(operator_strings)):
    sel_freq_op = sel_freq[operator_strings[i]]
    plt.bar(nfe_array, sel_freq_op, bottom=sel_freq_sum)
    sel_freq_sum = np.add(sel_freq_sum, sel_freq_op)
plt.xlabel('NFE')
plt.ylabel('Selection Frequency')
plt.legend(heur_labels)
plt.show()

    

                    
            