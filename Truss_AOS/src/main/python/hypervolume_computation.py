# -*- coding: utf-8 -*-
"""
Hypervolume computation post GA run in Java

@author: roshan94
"""
from PyGMO.util import *
import csv
import numpy as np

#### Read data from the appropriate csv file

# set to: true - to read results for fibre stiffness model run
#         false - to read results for truss model run
fibre_stiffness = True

# set to: true - if Epsilon MOEA was used
#         false - if AOS MOEA was used
eps_moea = True

filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\'
if eps_moea:
    fileloc = 'Epsilon MOEA Runs\\'
    if fibre_stiffness:
        filename = 'Fibre Stiffness code run results\\EpsilonMOEA_emoea0_fibrestiffness_fullpop.csv'
    else: 
        filename = 'Truss code run results\\EpsilonMOEA_emoea0_trussstiffness_fullpop.csv'
else:
    fileloc = 'AOS MOEA Runs\\'
    if fibre_stiffness:
        filename = 'Fibre Stiffness code run results\\EpsilonMOEA_emoea0_fibrestiffness_fullpop.csv'
    else: 
        filename = 'Truss code run results\\EpsilonMOEA_emoea0_trussstiffness_fullpop.csv' 

full_filepath = filepath + fileloc + filename

with open(full_filepath,newline='') as csvfile:
    data = [row for row in csv.reader(csvfile)]
    designs = ["" for x in range(len(data)-1)]
    num_func_evals = np.zeros((len(data)-1,1))
    pen_obj1 = np.zeros((len(data)-1,1))
    pen_obj2 = np.zeros((len(data)-1,1))
    feas_scores = np.zeros((len(data)-1,1))
    stab_scores = np.zeros((len(data)-1,1))
    for x in range(len(data)-1):
        designs[x] = data[x+1][0]
        num_func_evals[x] = int(data[x+1][1])
        pen_obj1[x] = float(data[x+1][2])
        pen_obj2[x] = float(data[x+1][3])
        feas_scores[x] = float(data[x+1][4])
        stab_scores[x] = float(data[x+1][5])
        
n_des = len(designs)
des_array = np.zeros((n_des,36))
for x in range(n_des):
    current_des = designs[x]
    for y in range(36):
        des_array[x][y] = int(current_des[y])

#### Compute true objectives
obj1 = np.zeros(n_des)
obj2 = np.zeros(n_des)
pen_fac = 1
if fibre_stiffness:
    pen_fac = 2

for x in range(n_des):
    pen = (np.log10(np.absolute(feas_scores[x])) + np.log10(np.absolute(stab_scores[x])))/2
    obj1[x] = 15*(pen_obj1[x] + pen_fac*pen)
    obj2[x] = -85000*(pen_obj2[x] + pen_fac*pen)

#### Compute and plot hypervolume values


