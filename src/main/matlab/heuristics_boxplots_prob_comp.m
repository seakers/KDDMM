%% Feasibility, Stability and Orientation boxplots for Interior Penalty - Feas case
clear
close all
clc

%% Program Operation
%%%% Read csv files
fibre_model_used = false;

num_runs = 30; % change based on run 

% feas_bools = [int_pen, AOS, biased_init, ACH] boolean array
% stab_bools = [int_pen, AOS, biased_init, ACH] boolean array
% orient_bools = [int_pen, AOS, biased_init, ACH] boolean array

% Case - int_pen - Feas
case_feas_bools = [true, false, false, false];
case_stab_bools = [false, false, false, false];
case_orient_bools = [false, false, false, false];

prob2_used = false; % boolean whether optimization problem 2 data is to be read

[feas_struct, stab_struct, orient_struct] = get_heuristic_scores_oa(fibre_model_used, case_feas_bools, case_stab_bools, case_orient_bools, prob2_used, num_runs);

labels = cell(num_runs,1);
for i = 1:num_runs
    labels{i,1} = num2str(i-1);
end

% heuristic box plots for case 
[feas_array_all, feas_array_mean, feas_groups] = create_boxplot_arrays(feas_struct, num_runs);
figure
boxplot(feas_array_all,feas_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Feasiblity Score')
title('Feasibliity Score Boxplot for each run')

figure 
boxplot(feas_array_mean);
%boxplot(feas_array_all);
ylabel('Feasiblity Score')
title('Feasibliity Score Boxplot over all runs')

[stab_array_all, stab_array_mean, stab_groups] = create_boxplot_arrays(stab_struct, num_runs);
figure
boxplot(stab_array_all,stab_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Stability Score')
title('Stability Score Boxplot for each run')

figure 
boxplot(stab_array_mean);
%boxplot(stab_array_all);
ylabel('Stability Score')
title('Stability Score Boxplot over all runs')

[orient_array_all, orient_array_mean, orient_groups] = create_boxplot_arrays(orient_struct, num_runs);
figure
boxplot(orient_array_all,orient_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Orientation Score')
title('Orientation Score Boxplot for each run')

figure 
boxplot(orient_array_mean);
%boxplot(orient_array_all);
ylabel('Orientation Score')
title('Orientation Score Boxplot overe all runs')

%% Functions
function [feas_struct_all, stab_struct_all, orient_struct_all] = get_heuristic_scores_oa(fib_stiff, feas_bools, stab_bools, orient_bools, prob2_bool, n_runs)
    feas_struct_all = struct;
    stab_struct_all = struct;
    orient_struct_all = struct;
    
    filepath = "C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\";
    methods = ["Interior Penalty","AOS","Biased Initialization","Adaptive Constraint Handling"];
    
    if prob2_bool
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
        filename_model = "_fibre.csv";
    else
        filepath3 = "Truss Stiffness\\";
        filename_model = "_truss.csv";
    end
      
    %%%% read appropriate files 
    for i = 1:n_runs
        run_num = i-1;
        
        full_filepath = strcat(filepath,filepath_prob,filepath2,filepath_moea,filepath3,filename,num2str(run_num),filename2,filename_prob,filename_model);

        data_table = readtable(full_filepath,'Format','%s%f%f%f%f%f%f%f','HeaderLines',1);

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
    
        f_penalized = csv_data_array(:,1:2);
        feas_array = csv_data_array(:,5);
        stab_array = csv_data_array(:,6);
        orient_array = csv_data_array(:,7);
        
        pareto_bool = paretofront(f_penalized);
        designs_pareto = designs_array(pareto_bool==1);
        feas_pareto = feas_array(pareto_bool==1);
        stab_pareto = stab_array(pareto_bool==1); 
        orient_pareto = orient_array(pareto_bool==1);
        
        current_field = strcat('run_',num2str(i));
        
        feas_struct_all.(current_field) = feas_pareto;
        stab_struct_all.(current_field) = stab_pareto;
        orient_struct_all.(current_field) = orient_pareto;
        
    end
end

function [val_array, mean_val_array, bp_groups] = create_boxplot_arrays(val_structs, n_runs)
    val_array = [];
    bp_groups = [];
    mean_val_array = zeros(n_runs,1);
    group_count = 0;
    for i = 1:n_runs
        current_field = strcat('run_',num2str(i));
        current_val_array = val_structs.(current_field);
        val_array = [val_array;current_val_array];
        bp_groups = [bp_groups;group_count.*ones(size(current_val_array,1),1)];
        mean_val_array(i) = mean(current_val_array);
        group_count = group_count + 1;
    end
end

