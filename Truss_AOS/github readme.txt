Scripts to generate datasets for the Metamaterial problems for the paper: 
Suresh Kumar, Roshan, Srikar Srivatsa, Emilie Baker, Meredith Silberstein, and Daniel Selva. "Identifying and Leveraging Promising Design Heuristics for Multi-Objective Combinatorial Design Optimization." Journal of Mechanical Design 145, no. 12 (2023).

Important scripts:
JAVA:
ConstantRadiusMOEARun.java - Start and store results for multiple runs of either problem with different heuristic implementations
GenerateOperatorIndexDataset.java - Generate datasets for repair operators screening study for either problem
GenerateBiasedSamplingDataset.java - Generate datasets for biased sampling screening study for either problem

MATLAB:
metrics_study_biasing.m - Conduct soft constraints screening study (random sampling datasets for screening study are generated within the script)

PYTHON:
hv_truss_material_heurcomp.py - Compute hypervolumes and statistics for different cases (efficacy study results)
operator_index_computation.py - Compute HDIs for repair operators 
biased_sampling_index_computation.py - Compute HDI for Orientation biased sampling function
