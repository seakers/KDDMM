# -*- coding: utf-8 -*-
"""
Test for combined constraint satisfaction and pareto dominance
function

@author: roshan94
"""
import numpy as np

#### Assumes that all objectives must be minimized
def compute_pareto_front(population_objs, population_constr_aggr):
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

#### Testing 

pop_objs = [[1,7],[3,8],[3,6],[4,6],[2,6],[3,5],[5,5],[4,4],[5,3],[4,8]]
#pop_constrs_aggr = list(np.zeros(len(pop_objs))) # Sum of absolute constraint violations
pop_constrs_aggr = [0.5,0.5,0,0,0.75,0.6,0.3,0.25,0,0.8] # Sum of absolute constraint violations

# pop_objs = [[1,6]]
# pop_constrs_aggr = [5]

pareto_objs = compute_pareto_front(pop_objs, pop_constrs_aggr)
print(pareto_objs)