# -*- coding: utf-8 -*-
"""
Plot credits and quality values for each heuristic for a single run of the 
metamaterial problems

@author: roshan94
"""
import numpy as np
import csv
import matplotlib.pyplot as plt

artery_problem = False
model = 1 # 0 -> Fibre model, 1 -> Truss model, 2 -> Beam model

### Read the credit and quality csv files
file_loc = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\ahs truss good\\'

ahs_heur_bools = [True, True, True, True] # [part_coll, nodal_prop, orient, inters]

run_num = 1 # run number of results to read

heur_abbrvs_list = ['p','n','o','i']

filename_prob = '_prob2'
if artery_problem:
    filename_prob = '_artery'
    
filename_model = '_fibre'
if model == 1:
    filename_model = '_truss'
elif model == 2:
    filename_model = '_beam'
    
heurs = ''
for i in range(len(ahs_heur_bools)):
    if ahs_heur_bools[i]:
        heurs = heurs + heur_abbrvs_list[i]

cred_filename = file_loc + 'emoea_' + str(run_num) + heurs + 'con6_' + filename_prob + filename_model + '_credit.csv'
qual_filename = file_loc + 'emoea_' + str(run_num) + heurs + 'con6_' + filename_prob + filename_model + '_qual.csv'

heurs_strings = ['PartialCollapsibilityViolation','NodalPropertiesViolation','OrientationViolation','IntersectionViolation']

cred_vals = {}
qual_vals = {}

with open(cred_filename,newline='') as csvfile:
    file_reader = csv.reader(csvfile)
    
    ind_row = 0
    for row in file_reader:
        if ((ind_row > 0) and (row[0] == 'iteration')):
            continue
        else:
            if (ind_row == 0):
                nfe_array_all = row[1:]
            if row[0] in heurs_strings:
                heur_ind = heurs_strings.index(row[0])
                heur = heurs_strings[heur_ind]
                cred_vals[heur] = row[1:]
        ind_row = ind_row + 1
    
with open(qual_filename,newline='') as csvfile:
    file_reader = csv.reader(csvfile)
    
    ind_row = 0
    for row in file_reader:
        if ((ind_row > 0) and (row[0] == 'iteration')):
            continue
        else:
            if row[0] in heurs_strings:
                heur_ind = heurs_strings.index(row[0])
                heur = heurs_strings[heur_ind]
                qual_vals[heur] = row[1:]
        ind_row = ind_row + 1

### Average credits and quality values over each NFE val 
nfe_array_inds = {} # dict of the form nfe_val:[start_ind, end_ind] in nfe_array_all

current_ind = 0
current_nfe = nfe_array_all[current_ind]

while (True):
    max_ind_current_nfe = max(index for index, item in enumerate(nfe_array_all) if item == current_nfe)
    nfe_array_inds[current_nfe] = [current_ind, max_ind_current_nfe]
    
    if (max_ind_current_nfe == (len(nfe_array_all)-1)):
        break
    
    current_ind = max_ind_current_nfe + 1
    current_nfe = nfe_array_all[current_ind]
    
nfe_array = []
creds_arrays = {}
quals_arrays = {}

nfe_array_populated = False

for i in range(len(heurs_strings)):
    cred_vals_heur = np.array(cred_vals[heurs_strings[i]]).astype(float)
    qual_vals_heur = np.array(qual_vals[heurs_strings[i]]).astype(float)
    
    cred_avg = []
    qual_avg = []
    
    for nfe in nfe_array_inds:
        if (not nfe_array_populated):
            nfe_array.append(int(nfe))
        nfe_inds = nfe_array_inds[nfe]
        
        cred_avg.append(np.mean(cred_vals_heur[nfe_inds[0]:nfe_inds[1]]))
        qual_avg.append(np.mean(qual_vals_heur[nfe_inds[0]:nfe_inds[1]]))
        
    creds_arrays[heurs_strings[i]] = cred_avg
    quals_arrays[heurs_strings[i]] = qual_avg
    nfe_array_populated = True
    
### Plot credits and quality values

# Credit plots
fig1 = plt.figure()
for i in range(len(heurs_strings)):
    cred_array = creds_arrays[heurs_strings[i]]
    plt.plot(nfe_array, cred_array, label=heurs_strings[i])
plt.xlabel('NFE')
plt.ylabel('Credit Value')
plt.legend(loc="upper right")
plt.show()

# Quality plots
fig2 = plt.figure()
for i in range(len(heurs_strings)):
    qual_array = quals_arrays[heurs_strings[i]]
    plt.plot(nfe_array, qual_array, label=heurs_strings[i])
plt.xlabel('NFE')
plt.ylabel('Quality Value')
plt.legend(loc="upper right")
plt.show()
    
    
    
        
        
    
                
            
    
    

