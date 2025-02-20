%% Check bounds of objective values
clear
close all
clc

%% Program Operation
%%%% Read csv files
fib_stiff = false;

n_runs = 30; % change based on run 

% feas_bools = [int_pen, AOS, biased_init, ACH] boolean array
% stab_bools = [int_pen, AOS, biased_init, ACH] boolean array
% orient_bools = [int_pen, AOS, biased_init, ACH] boolean array

feas_bools = [false, false, false, false];
stab_bools = [false, false, false, false];
orient_bools = [false, false, false, false];

prob2_used = false; % boolean whether optimization problem 2 data is to be read

filepath = "C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\";
methods = ["Interior Penalty","AOS","Biased Initialization","Adaptive Constraint Handling"];

if prob2_used
    filepath_prob = "OA run data - Optimization Problem 2\\";
    filename_prob = "_prob2";
else
    filepath_prob = "OA run data - Optimization Problem 1\\";
    filename_prob = "";
end
    
if (feas_bools(2) || stab_bools(2) || orient_bools(2))
    filename = "AOSMOEA_emoea_";
else
    filename = "EpsilonMOEA_emoea_";
end
filepath2 = '';
filename2 = '';
constr_count = 0;
for i = 1:4
    if (feas_bools(i) && stab_bools(i) && ~orient_bools(i))
        constrained = strcat(methods(i)," - Feasibility and Stability\\");
        filename_snippet = strcat("fscon",num2str(i-1),"_");
    elseif (feas_bools(i) && ~stab_bools(i) && ~orient_bools(i))
        constrained = strcat(methods(i)," - Feasibility\\");
        filename_snippet = strcat("fcon",num2str(i-1),"_");
    elseif (~feas_bools(i) && stab_bools(i) && ~orient_bools(i))
        constrained = strcat(methods(i)," - Stability\\");
        filename_snippet = strcat("scon",num2str(i-1),"_");
    elseif (~feas_bools(i) && ~stab_bools(i) && orient_bools(i))
        constrained = strcat(methods(i)," - Orientation\\");
        filename_snippet = strcat("ocon",num2str(i-1),"_");
    elseif (feas_bools(i) && ~stab_bools(i) && orient_bools(i))
        constrained = strcat(methods(i)," - Feasibility and Orientation\\");
        filename_snippet = strcat("focon",num2str(i-1),"_");
    elseif (~feas_bools(i) && stab_bools(i) && orient_bools(i))
        constrained = strcat(methods(i)," - Stability and Orientation\\");
        filename_snippet = strcat("socon",num2str(i-1),"_");
    elseif (feas_bools(i) && stab_bools(i) && orient_bools(i))
        constrained = strcat(methods(i)," - Feasiblity Stability and Orientation\\");
        filename_snippet = strcat("fsocon",num2str(i-1),"_");
    else
        constrained = '';
        filename_snippet = '';
        constr_count = constr_count + 1;
    end
    filepath2 = strcat(filepath2, constrained);
    filename2 = strcat(filename2, filename_snippet);
end
    
filepath_moea = '';
if (constr_count == 4)
    filepath_moea = "Epsilon MOEA\\";
end

if fib_stiff
    filepath3 = "Fibre Stiffness\\";
    filename_model = "_fibre_fullpop.csv";
else
    filepath3 = "Truss Stiffness\\";
    filename_model = "_truss_fullpop.csv";
end

%%%% read appropriate files 
f_pen_bounds = zeros(2,2,n_runs);
f_true_bounds = zeros(2,2,n_runs);
for i = 1:n_runs
    run_num = i-1;
    
    full_filepath = strcat(filepath,filepath_prob,filepath2,filepath_moea,filepath3,filename,num2str(run_num),filename2,filename_prob,filename_model);
    
    data_table = readtable(full_filepath,'Format','%s%f%f%f%f%f%f%f%f','HeaderLines',1,'ReadVariableNames',false);
    
    %%%% store retrieved data into different variables
    %%%% csv_data includes: [Pen. Obj. 1, Pen.Obj. 2, Feasibility Score,
    %%%% Stablity Score]
    pop_size =  size(data_table,1);
    csv_data = zeros(pop_size,7);
    designs = strings(pop_size);
    csv_data = data_table(:,2:end);
    designs = data_table(:,1);
    
    csv_data_array = table2array(csv_data);
    designs_array = table2array(designs);
    
    f_penalized = csv_data_array(:,2:3);
    f_true = csv_data_array(:,4:5);

    f_pen_bounds(:,:,i) = [min(f_penalized(:,1)),max(f_penalized(:,1)); min(f_penalized(:,2)), max(f_penalized(:,2))];
    f_true_bounds(:,:,i) = [min(f_true(:,1)),max(f_true(:,1)); min(f_true(:,2)), max(f_true(:,2))];
    
end

%% Find bounds on objectives over all runs
f1_pen_bounds_allruns = zeros(n_runs,2);
f2_pen_bounds_allruns = zeros(n_runs,2);
for i = 1:n_runs
    f1_pen_bounds_allruns(i,:) = f_pen_bounds(1,:,i);
    f2_pen_bounds_allruns(i,:) = f_pen_bounds(2,:,i);
end 

f1_true_bounds_allruns = zeros(n_runs,2);
f2_true_bounds_allruns = zeros(n_runs,2);
for i = 1:n_runs
    f1_true_bounds_allruns(i,:) = f_true_bounds(1,:,i);
    f2_true_bounds_allruns(i,:) = f_true_bounds(2,:,i);
end 

f1_pen_bounds_aggregate = [min(f1_pen_bounds_allruns(:,1)), max(f1_pen_bounds_allruns(:,2))];
f2_pen_bounds_aggregate = [min(f2_pen_bounds_allruns(:,1)), max(f2_pen_bounds_allruns(:,2))];

f1_true_bounds_aggregate = [min(f1_true_bounds_allruns(:,1)), max(f1_true_bounds_allruns(:,2))];
f2_true_bounds_aggregate = [min(f2_true_bounds_allruns(:,1)), max(f2_true_bounds_allruns(:,2))];


