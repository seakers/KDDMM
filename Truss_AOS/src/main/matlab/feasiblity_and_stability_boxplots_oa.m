%% Aggregated Feasibility and Stability comparison (Orthogonal Array runs)
clear
close all
clc

%% Program Operation
%%%% Read csv files
fibre_model_used = false;

% feas_bools = [int_pen, AOS, biased_init, ACH] boolean array
% stab_bools = [int_pen, AOS, biased_init, ACH] boolean array

% Case 1 - Epsilon MOEA
case1_feas_bools = [false, false, false, false];
case1_stab_bools = [false, false, false, false];

% Case 2 - bias_init - Feas
case2_feas_bools = [false, false, true, false];
case2_stab_bools = [false, false, false, false];

% Case 3 - int_pen - Feas
case3_feas_bools = [true, false, false, false];
case3_stab_bools = [false, false, false, false];

% Case 4 - bias_init - Feas, AOS - Stab
case4_feas_bools = [false, false, true, false];
case4_stab_bools = [false, true, false, false];

% Case 5 - int_pen - Feas, ACH - Stab
case5_feas_bools = [true, false, false, false];
case5_stab_bools = [false, false, false, true];

% Case 6 - biased_int & int_pen - Feas
case6_feas_bools = [true, false, true, false];
case6_stab_bools = [false, false, false, false];

% Case 7 - int_pen - Feas & Stab
case7_feas_bools = [true, false, false, false];
case7_stab_bools = [true, false, false, false];

% Case 8 - biased_init & AOS - Feas
case8_feas_bools = [false, true, true, false];
case8_stab_bools = [false, false, false, false];

% Case 9 - int_pen - Feas, AOS - Stab
case9_feas_bools = [true, false, false, false];
case9_stab_bools = [false, true, false, false];

num_runs = 30; % change based on run 
pop_size = 100; % change based on run

[feas_struct_case1, stab_struct_case1] = get_feas_and_stab_scores_oa(fibre_model_used, case1_feas_bools, case1_stab_bools, num_runs, pop_size);
[feas_struct_case2, stab_struct_case2] = get_feas_and_stab_scores_oa(fibre_model_used, case2_feas_bools, case2_stab_bools, num_runs, pop_size);
[feas_struct_case3, stab_struct_case3] = get_feas_and_stab_scores_oa(fibre_model_used, case3_feas_bools, case3_stab_bools, num_runs, pop_size);
[feas_struct_case4, stab_struct_case4] = get_feas_and_stab_scores_oa(fibre_model_used, case4_feas_bools, case4_stab_bools, num_runs, pop_size);
[feas_struct_case5, stab_struct_case5] = get_feas_and_stab_scores_oa(fibre_model_used, case5_feas_bools, case5_stab_bools, num_runs, pop_size);
[feas_struct_case6, stab_struct_case6] = get_feas_and_stab_scores_oa(fibre_model_used, case6_feas_bools, case6_stab_bools, num_runs, pop_size);
[feas_struct_case7, stab_struct_case7] = get_feas_and_stab_scores_oa(fibre_model_used, case7_feas_bools, case7_stab_bools, num_runs, pop_size);
[feas_struct_case8, stab_struct_case8] = get_feas_and_stab_scores_oa(fibre_model_used, case8_feas_bools, case8_stab_bools, num_runs, pop_size);
[feas_struct_case9, stab_struct_case9] = get_feas_and_stab_scores_oa(fibre_model_used, case9_feas_bools, case9_stab_bools, num_runs, pop_size);

%% Create boxplots for individual cases
labels = cell(num_runs,1);
for i = 1:num_runs
    labels{i,1} = num2str(i-1);
end

% feasibility and stability box plots for case 1
[feas_array_all_case1, feas_array_mean_case1, feas_case1_groups] = create_boxplot_arrays(feas_struct_case1, num_runs);
figure
boxplot(feas_array_all_case1,feas_case1_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Feasiblity Score')
title('Feasibliity Score Boxplot for Case 1')

[stab_array_all_case1, stab_array_mean_case1, stab_case1_groups] = create_boxplot_arrays(stab_struct_case1, num_runs);
figure
boxplot(stab_array_all_case1,stab_case1_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Stability Score')
title('Stability Score Boxplot for Case 1')

% feasibility and stability box plots for case 2
[feas_array_all_case2, feas_array_mean_case2, feas_case2_groups] = create_boxplot_arrays(feas_struct_case2, num_runs);
figure
boxplot(feas_array_all_case2,feas_case2_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Feasiblity Score')
title('Feasibliity Score Boxplot for Case 2')

[stab_array_all_case2, stab_array_mean_case2, stab_case2_groups] = create_boxplot_arrays(stab_struct_case2, num_runs);
figure
boxplot(stab_array_all_case2,stab_case2_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Stability Score')
title('Stability Score Boxplot for Case 2')

% feasibility and stability box plots for case 3
[feas_array_all_case3, feas_array_mean_case3, feas_case3_groups] = create_boxplot_arrays(feas_struct_case3, num_runs);
figure
boxplot(feas_array_all_case3,feas_case3_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Feasiblity Score')
title('Feasibliity Score Boxplot for Case 3')

[stab_array_all_case3, stab_array_mean_case3, stab_case3_groups] = create_boxplot_arrays(stab_struct_case3, num_runs);
figure
boxplot(stab_array_all_case3,stab_case3_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Stability Score')
title('Stability Score Boxplot for Case 3')

% feasibility and stability box plots for case 4
[feas_array_all_case4, feas_array_mean_case4, feas_case4_groups] = create_boxplot_arrays(feas_struct_case4, num_runs);
figure
boxplot(feas_array_all_case4,feas_case4_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Feasiblity Score')
title('Feasibliity Score Boxplot for Case 4')

[stab_array_all_case4, stab_array_mean_case4, stab_case4_groups] = create_boxplot_arrays(stab_struct_case4, num_runs);
figure
boxplot(stab_array_all_case4,stab_case4_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Stability Score')
title('Stability Score Boxplot for Case 4')

% feasibility and stability box plots for case 5
[feas_array_all_case5, feas_array_mean_case5, feas_case5_groups] = create_boxplot_arrays(feas_struct_case5, num_runs);
figure
boxplot(feas_array_all_case5,feas_case5_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Feasiblity Score')
title('Feasibliity Score Boxplot for Case 5')

[stab_array_all_case5, stab_array_mean_case5, stab_case5_groups] = create_boxplot_arrays(stab_struct_case5, num_runs);
figure
boxplot(stab_array_all_case5,stab_case5_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Stability Score')
title('Stability Score Boxplot for Case 5')

% feasibility and stability box plots for case 6
[feas_array_all_case6, feas_array_mean_case6, feas_case6_groups] = create_boxplot_arrays(feas_struct_case6, num_runs);
figure
boxplot(feas_array_all_case6,feas_case6_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Feasiblity Score')
title('Feasibliity Score Boxplot for Case 6')

[stab_array_all_case6, stab_array_mean_case6, stab_case6_groups] = create_boxplot_arrays(stab_struct_case6, num_runs);
figure
boxplot(stab_array_all_case6,stab_case6_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Stability Score')
title('Stability Score Boxplot for Case 6')

% feasibility and stability box plots for case 7
[feas_array_all_case7, feas_array_mean_case7, feas_case7_groups] = create_boxplot_arrays(feas_struct_case7, num_runs);
figure
boxplot(feas_array_all_case7,feas_case7_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Feasiblity Score')
title('Feasibliity Score Boxplot for Case 7')

[stab_array_all_case7, stab_array_mean_case7, stab_case7_groups] = create_boxplot_arrays(stab_struct_case7, num_runs);
figure
boxplot(stab_array_all_case7,stab_case7_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Stability Score')
title('Stability Score Boxplot for Case 7')

% feasibility and stability box plots for case 8
[feas_array_all_case8, feas_array_mean_case8, feas_case8_groups] = create_boxplot_arrays(feas_struct_case8, num_runs);
figure
boxplot(feas_array_all_case8,feas_case8_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Feasiblity Score')
title('Feasibliity Score Boxplot for Case 8')

[stab_array_all_case8, stab_array_mean_case8, stab_case8_groups] = create_boxplot_arrays(stab_struct_case8, num_runs);
figure
boxplot(stab_array_all_case8,stab_case8_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Stability Score')
title('Stability Score Boxplot for Case 8')

% feasibility and stability box plots for case 9
[feas_array_all_case9, feas_array_mean_case9, feas_case9_groups] = create_boxplot_arrays(feas_struct_case9, num_runs);
figure
boxplot(feas_array_all_case9,feas_case9_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Feasiblity Score')
title('Feasibliity Score Boxplot for Case 9')

[stab_array_all_case9, stab_array_mean_case9, stab_case9_groups] = create_boxplot_arrays(stab_struct_case9, num_runs);
figure
boxplot(stab_array_all_case9,stab_case9_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Stability Score')
title('Stability Score Boxplot for Case 9')

%% Create boxplots camparing different cases
case_labels = {'Case 1','Case 2','Case 3','Case 4','Case 5','Case 6','Case 7','Case 8','Case 9'};
mean_bp_groups = [zeros(1,num_runs),ones(1,num_runs),2.*ones(1,num_runs),3.*ones(1,num_runs),4.*ones(1,num_runs),5.*ones(1,num_runs),6.*ones(1,num_runs),7.*ones(1,num_runs),8.*ones(1,num_runs)];

% Plotting boxplots
feas_array = [feas_array_mean_case1',feas_array_mean_case2',feas_array_mean_case3',feas_array_mean_case4',feas_array_mean_case5',feas_array_mean_case6',feas_array_mean_case7',feas_array_mean_case8',feas_array_mean_case9'];
figure 
boxplot(feas_array,mean_bp_groups,'Labels',case_labels);
ylabel('Feasiblity Score')
title('Feasibliity Score Comparison Boxplot')

stab_array = [stab_array_mean_case1',stab_array_mean_case2',stab_array_mean_case3',stab_array_mean_case4',stab_array_mean_case5',stab_array_mean_case6',stab_array_mean_case7',stab_array_mean_case8',stab_array_mean_case9'];
figure 
boxplot(stab_array,mean_bp_groups,'Labels',case_labels);
ylabel('Stability Score')
title('Stability Score Comparison Boxplot')

%% Functions
function [feas_struct_all, stab_struct_all] = get_feas_and_stab_scores_oa(fib_stiff, feas_bools, stab_bools, n_runs, n_pop)
    feas_struct_all = struct;
    stab_struct_all = struct;
    
    filepath = "C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\";
    methods = ["Interior Penalty","AOS","Biased Initialization","Adaptive Constraint Handling"];
    
    if (feas_bools(2) || stab_bools(2))
        filename = "AOSMOEA_emoea_";
    else
        filename = "EpsilonMOEA_emoea_";
    end
    filepath2 = '';
    filename2 = '';
    constr_count = 0;
    for i = 1:4
        if (feas_bools(i) && stab_bools(i))
            constrained = strcat(methods(i)," - Feasibility and Stability\\");
            filename_snippet = strcat("fscon",num2str(i-1),"_");
        elseif (feas_bools(i) && ~stab_bools(i))
            constrained = strcat(methods(i)," - Feasibility\\");
            filename_snippet = strcat("fcon",num2str(i-1),"_");
        elseif (~feas_bools(i) && stab_bools(i))
            constrained = strcat(methods(i)," - Stability\\");
            filename_snippet = strcat("scon",num2str(i-1),"_");
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
        filename2 = strcat(filename2,"_fibre.csv");
    else
        filepath3 = "Truss Stiffness\\";
        filename2 = strcat(filename2,"_truss.csv");
    end
      
    %%%% read appropriate files 
    for i = 1:n_runs
        run_num = i-1;
        
        full_filepath = strcat(filepath,filepath2,filepath_moea,filepath3,filename,num2str(run_num),filename2);

        data_table = readtable(full_filepath,'Format','%s%f%f%f%f%f%f','HeaderLines',1);

        %%%% store retrieved data into different variables
        %%%% csv_data includes: [Pen. Obj. 1, Pen.Obj. 2, Feasibility Score,
        %%%% Stablity Score]
        pop_size =  size(data_table,1);
        csv_data = zeros(pop_size,6);
        designs = strings(pop_size);
        csv_data = data_table(:,2:end);
        designs = data_table(:,1);
    
        csv_data_array = table2array(csv_data);
        designs_array = table2array(designs);
    
        f_penalized = csv_data_array(:,1:2);
        feas_array = csv_data_array(:,5);
        stab_array = csv_data_array(:,6);
        
        pareto_bool = paretofront(f_penalized);
        designs_pareto = designs_array(pareto_bool==1);
        feas_pareto = feas_array(pareto_bool==1);
        stab_pareto = stab_array(pareto_bool==1); 
        
        current_field = strcat('run_',num2str(i));
        
        feas_struct_all.(current_field) = feas_pareto;
        stab_struct_all.(current_field) = stab_pareto;
        
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