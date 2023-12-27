# KDDMM
Design optimization tool for meta-materials that facilitates the fusion of expert knowledge (physics-based models, heuristics) and data-driven approaches (surrogate models). 

Scripts to generate datasets for the Metamaterial problems for the paper: 
Suresh Kumar, Roshan, Srikar Srivatsa, Emilie Baker, Meredith Silberstein, and Daniel Selva. "Identifying and Leveraging Promising Design Heuristics for Multi-Objective Combinatorial Design Optimization." Journal of Mechanical Design 145, no. 12 (2023).

## Dependencies:
### JAVA (can be found in the Maven pom.xml file):
1. MOEAFramework: https://github.com/MOEAFramework/MOEAFramework
2. mopAOS: https://github.com/seakers/mopAOS/tree/heuristics (heuristics branch)
3. Adaptive Heuristic Selection: https://github.com/seakers/Adaptive-Heuristic-Selection
4. System Architecture Problems: https://github.com/seakers/SystemArchitectureProblems
5. Mathworks engine (for design evaluation): R2020a is used, change according to your Matlab version

### Python:
1. PyGMO (Python Parallel Global Multiobjective Optimizer): https://esa.github.io/pygmo/
2. Scipy
3. Matplotlib

## Important scripts:
### JAVA:
1. ConstantRadiusMOEARun.java - Start and store results for multiple runs of either problem with different heuristic implementations
2. GenerateOperatorIndexDataset.java - Generate datasets for repair operators screening study for either problem
3. GenerateBiasedSamplingDataset.java - Generate datasets for biased sampling screening study for either problem

### MATLAB:
1. metrics_study_biasing.m - Conduct soft constraints screening study (random sampling datasets for screening study are generated within the script)

### PYTHON:
1. hv_truss_material_heurcomp.py - Compute hypervolumes and statistics for different cases (efficacy study results)
2. operator_index_computation.py - Compute HDIs for repair operators 
3. biased_sampling_index_computation.py - Compute HDI for Orientation biased sampling function
