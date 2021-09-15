%% Aggregated Heuristics and Constraints comparison
clear
close all
clc

%% Program Parameters
%%%% Read csv files
fibre_model_used = false;
sidenum = 3;

num_runs = 30; % change based on run 
pop_size = 100; % change based on run

% partcoll_bools = [int_pen, AOS, biased_init, ACH] boolean array
% nodalprop_bools = [int_pen, AOS, biased_init, ACH] boolean array
% orient_bools = [int_pen, AOS, biased_init, ACH] boolean array

%% Create boxplots comparing different cases (eps-moea vs AOS - Orientation)(Constant Radii problem cases)
constrad_read = true;
% Case 1 - Epsilon MOEA
repeat3x3_case1 = true;
case1_partcoll_bools = [false, false, false, false];
case1_nodalprop_bools = [false, false, false, false];
case1_orient_bools = [false, false, false, false];
case1_feas_bools = [false, false];

% Case 2 - AOS - Orientation 
repeat3x3_case2 = true;
case2_partcoll_bools = [false, false, false, false];
case2_nodalprop_bools = [false, false, false, false];
case2_orient_bools = [false, true, false, false];
case2_feas_bools = [true, true];

[feas_struct_case1, conn_struct_case1, stiffrat_struct_case1, partcoll_struct_case1, nodalprop_struct_case1, orient_struct_case1] = get_heuristic_scores_oa(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_feas_bools, constrad_read, repeat3x3_case1, sidenum, num_runs, pop_size);
[feas_struct_case2, conn_struct_case2, stiffrat_struct_case2, partcoll_struct_case2, nodalprop_struct_case2, orient_struct_case2] = get_heuristic_scores_oa(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_feas_bools, constrad_read, repeat3x3_case2, sidenum, num_runs, pop_size);

labels = cell(num_runs,1);
for i = 1:num_runs
    labels{i,1} = num2str(i-1);
end

% heuristic and constraint box plots for case 1
[feas_array_all_case1, feas_array_mean_case1, feas_case1_groups] = create_boxplot_arrays(feas_struct_case1, num_runs);
[conn_array_all_case1, conn_array_mean_case1, conn_case1_groups] = create_boxplot_arrays(conn_struct_case1, num_runs);
[stiffrat_array_all_case1, stiffrat_array_mean_case1, stiffrat_case1_groups] = create_boxplot_arrays(stiffrat_struct_case1, num_runs);
[partcoll_array_all_case1, partcoll_array_mean_case1, partcoll_case1_groups] = create_boxplot_arrays(partcoll_struct_case1, num_runs);
[nodalprop_array_all_case1, nodalprop_array_mean_case1, nodalprop_case1_groups] = create_boxplot_arrays(nodalprop_struct_case1, num_runs);
[orient_array_all_case1, orient_array_mean_case1, orient_case1_groups] = create_boxplot_arrays(orient_struct_case1, num_runs);

% heuristic and constraint box plots for case 2
[feas_array_all_case2, feas_array_mean_case2, feas_case2_groups] = create_boxplot_arrays(feas_struct_case2, num_runs);
[conn_array_all_case2, conn_array_mean_case2, conn_case2_groups] = create_boxplot_arrays(conn_struct_case2, num_runs);
[stiffrat_array_all_case2, stiffrat_array_mean_case2, stiffrat_case2_groups] = create_boxplot_arrays(stiffrat_struct_case2, num_runs);
[partcoll_array_all_case2, partcoll_array_mean_case2, partcoll_case2_groups] = create_boxplot_arrays(partcoll_struct_case2, num_runs);
[nodalprop_array_all_case2, nodalprop_array_mean_case2, nodalprop_case2_groups] = create_boxplot_arrays(nodalprop_struct_case2, num_runs);
[orient_array_all_case2, orient_array_mean_case2, orient_case2_groups] = create_boxplot_arrays(orient_struct_case2, num_runs);

case_labels = {'Without Heuristics','With Heuristics'};
mean_bp_groups = [zeros(1,num_runs),ones(1,num_runs)];

% Plotting case comparison boxplots
feas_array = [feas_array_mean_case1',feas_array_mean_case2'];
figure 
boxplot(feas_array,mean_bp_groups,'Labels',case_labels);
%feas_array = [feas_array_all_case1',feas_array_all_case2'];
%bp_groups_feas = [zeros(1,size(feas_array_all_case1,1)),ones(1,size(feas_array_all_case2,1))];
%boxplot(feas_array,bp_groups_feas,'Labels',case_labels);
ylabel('Feasiblity Score')
%title('Feasibliity Score Comparison Boxplot')

conn_array = [conn_array_mean_case1',conn_array_mean_case2'];
figure 
boxplot(conn_array,mean_bp_groups,'Labels',case_labels);
%conn_array = [conn_array_all_case1',conn_array_all_case2'];
%bp_groups_conn = [zeros(1,size(conn_array_all_case1,1)),ones(1,size(conn_array_all_case2,1))];
%boxplot(conn_array,bp_groups_conn,'Labels',case_labels);
ylabel('Connectivity Score')
%title('Connectivity Score Comparison Boxplot')

stiffrat_array = [stiffrat_array_mean_case1',stiffrat_array_mean_case2'];
figure 
boxplot(stiffrat_array,mean_bp_groups,'Labels',case_labels);
%stiffrat_array = [stiffrat_array_all_case1',stiffrat_array_all_case2'];
%bp_groups_stiffrat = [zeros(1,size(stiffrat_array_all_case1,1)),ones(1,size(stiffrat_array_all_case2,1))];
%boxplot(stiffrat_array,bp_groups_stiffrat,'Labels',case_labels);
ylabel('Stiffness Ratio Constraint')
%title('Stiffness Ratio Constraint Comparison Boxplot')

partcoll_array = [partcoll_array_mean_case1',partcoll_array_mean_case2'];
figure 
boxplot(partcoll_array,mean_bp_groups,'Labels',case_labels);
%partcoll_array = [partcoll_array_all_case1',partcoll_array_all_case2'];
%bp_groups_partcoll = [zeros(1,size(partcoll_array_all_case1,1)),ones(1,size(partcoll_array_all_case2,1))];
%boxplot(partcoll_array,bp_groups_partcoll,'Labels',case_labels);
ylabel('Partial Collapsibility Score')
%title('Partial Collapsibility Score Comparison Boxplot')

nodalprop_array = [nodalprop_array_mean_case1',nodalprop_array_mean_case2'];
figure 
boxplot(nodalprop_array,mean_bp_groups,'Labels',case_labels);
%nodalprop_array = [nodalprop_array_all_case1',nodalprop_array_all_case2'];
%bp_groups_nodalprop = [zeros(1,size(nodalprop_array_all_case1,1)),ones(1,size(nodalprop_array_all_case2,1))];
%boxplot(nodalprop_array,bp_groups_nodalprop,'Labels',case_labels);
ylabel('Nodal Properties Score')
%title('Nodal Properties Score Comparison Boxplot')

orient_array = [orient_array_mean_case1',orient_array_mean_case2'];
figure 
boxplot(orient_array,mean_bp_groups,'Labels',case_labels);
%orient_array = [orient_array_all_case1',orient_array_all_case2'];
%bp_groups_orient = [zeros(1,size(orient_array_all_case1,1)),ones(1,size(orient_array_all_case2,1))];
%boxplot(orient_array,bp_groups_orient,'Labels',case_labels);
ylabel('Orientation Score')
%title('Orientation Score Comparison Boxplot')

%% Create boxplots comparing different cases (eps-moea vs AOS - Orientation)(Variable Radii problem cases)
constrad_read = false;
% Case 1 - Epsilon MOEA
repeat3x3_case1 = true;
case1_partcoll_bools = [false, false, false, false];
case1_nodalprop_bools = [false, false, false, false];
case1_orient_bools = [false, false, false, false];
case1_feas_bools = [false, false];

% Case 2 - AOS - Orientation
repeat3x3_case2 = true;
case2_partcoll_bools = [false, false, false, false];
case2_nodalprop_bools = [false, false, false, false];
case2_orient_bools = [false, true, false, false];
case2_feas_bools = [true, true];

[feas_struct_case1, conn_struct_case1, stiffrat_struct_case1, partcoll_struct_case1, nodalprop_struct_case1, orient_struct_case1] = get_heuristic_scores_oa(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_feas_bools, constrad_read, repeat3x3_case1, sidenum, num_runs, pop_size);
[feas_struct_case2, conn_struct_case2, stiffrat_struct_case2, partcoll_struct_case2, nodalprop_struct_case2, orient_struct_case2] = get_heuristic_scores_oa(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_feas_bools, constrad_read, repeat3x3_case2, sidenum, num_runs, pop_size);

labels = cell(num_runs,1);
for i = 1:num_runs
    labels{i,1} = num2str(i-1);
end

% heuristic and constraint box plots for case 1
[feas_array_all_case1, feas_array_mean_case1, feas_case1_groups] = create_boxplot_arrays(feas_struct_case1, num_runs);
[conn_array_all_case1, conn_array_mean_case1, conn_case1_groups] = create_boxplot_arrays(conn_struct_case1, num_runs);
[stiffrat_array_all_case1, stiffrat_array_mean_case1, stiffrat_case1_groups] = create_boxplot_arrays(stiffrat_struct_case1, num_runs);
[partcoll_array_all_case1, partcoll_array_mean_case1, partcoll_case1_groups] = create_boxplot_arrays(partcoll_struct_case1, num_runs);
[nodalprop_array_all_case1, nodalprop_array_mean_case1, nodalprop_case1_groups] = create_boxplot_arrays(nodalprop_struct_case1, num_runs);
[orient_array_all_case1, orient_array_mean_case1, orient_case1_groups] = create_boxplot_arrays(orient_struct_case1, num_runs);

% heuristic and constraint box plots for case 2
[feas_array_all_case2, feas_array_mean_case2, feas_case2_groups] = create_boxplot_arrays(feas_struct_case2, num_runs);
[conn_array_all_case2, conn_array_mean_case2, conn_case2_groups] = create_boxplot_arrays(conn_struct_case2, num_runs);
[stiffrat_array_all_case2, stiffrat_array_mean_case2, stiffrat_case2_groups] = create_boxplot_arrays(stiffrat_struct_case2, num_runs);
[partcoll_array_all_case2, partcoll_array_mean_case2, partcoll_case2_groups] = create_boxplot_arrays(partcoll_struct_case2, num_runs);
[nodalprop_array_all_case2, nodalprop_array_mean_case2, nodalprop_case2_groups] = create_boxplot_arrays(nodalprop_struct_case2, num_runs);
[orient_array_all_case2, orient_array_mean_case2, orient_case2_groups] = create_boxplot_arrays(orient_struct_case2, num_runs);

case_labels = {'Without Heuristics','With Heuristics'};
mean_bp_groups = [zeros(1,num_runs),ones(1,num_runs)];

% Plotting case comparison boxplots
feas_array = [feas_array_mean_case1',feas_array_mean_case2'];
figure 
boxplot(feas_array,mean_bp_groups,'Labels',case_labels);
%feas_array = [feas_array_all_case1',feas_array_all_case2'];
%bp_groups_feas = [zeros(1,size(feas_array_all_case1,1)),ones(1,size(feas_array_all_case2,1))];
%boxplot(feas_array,bp_groups_feas,'Labels',case_labels);
ylabel('Feasiblity Score')
%title('Feasibliity Score Comparison Boxplot')

conn_array = [conn_array_mean_case1',conn_array_mean_case2'];
figure 
boxplot(conn_array,mean_bp_groups,'Labels',case_labels);
%conn_array = [conn_array_all_case1',conn_array_all_case2'];
%bp_groups_conn = [zeros(1,size(conn_array_all_case1,1)),ones(1,size(conn_array_all_case2,1))];
%boxplot(conn_array,bp_groups_conn,'Labels',case_labels);
ylabel('Connectivity Score')
%title('Connectivity Score Comparison Boxplot')

stiffrat_array = [stiffrat_array_mean_case1',stiffrat_array_mean_case2'];
figure 
boxplot(stiffrat_array,mean_bp_groups,'Labels',case_labels);
%stiffrat_array = [stiffrat_array_all_case1',stiffrat_array_all_case2'];
%bp_groups_stiffrat = [zeros(1,size(stiffrat_array_all_case1,1)),ones(1,size(stiffrat_array_all_case2,1))];
%boxplot(stiffrat_array,bp_groups_stiffrat,'Labels',case_labels);
ylabel('Stiffness Ratio Constraint')
%title('Stiffness Ratio Constraint Comparison Boxplot')

partcoll_array = [partcoll_array_mean_case1',partcoll_array_mean_case2'];
figure 
boxplot(partcoll_array,mean_bp_groups,'Labels',case_labels);
%partcoll_array = [partcoll_array_all_case1',partcoll_array_all_case2'];
%bp_groups_partcoll = [zeros(1,size(partcoll_array_all_case1,1)),ones(1,size(partcoll_array_all_case2,1))];
%boxplot(partcoll_array,bp_groups_partcoll,'Labels',case_labels);
ylabel('Partial Collapsibility Score')
%title('Partial Collapsibility Score Comparison Boxplot')

nodalprop_array = [nodalprop_array_mean_case1',nodalprop_array_mean_case2'];
figure 
boxplot(nodalprop_array,mean_bp_groups,'Labels',case_labels);
%nodalprop_array = [nodalprop_array_all_case1',nodalprop_array_all_case2'];
%bp_groups_nodalprop = [zeros(1,size(nodalprop_array_all_case1,1)),ones(1,size(nodalprop_array_all_case2,1))];
%boxplot(nodalprop_array,bp_groups_nodalprop,'Labels',case_labels);
ylabel('Nodal Properties Score')
%title('Nodal Properties Score Comparison Boxplot')

orient_array = [orient_array_mean_case1',orient_array_mean_case2'];
figure 
boxplot(orient_array,mean_bp_groups,'Labels',case_labels);
%orient_array = [orient_array_all_case1',orient_array_all_case2'];
%bp_groups_orient = [zeros(1,size(orient_array_all_case1,1)),ones(1,size(orient_array_all_case2,1))];
%boxplot(orient_array,bp_groups_orient,'Labels',case_labels);
ylabel('Orientation Score')
%title('Orientation Score Comparison Boxplot')

%% Functions
function [feas_struct_all, conn_struct_all, stiffrat_struct_all, partcoll_struct_all, nodalprop_struct_all, orient_struct_all] = get_heuristic_scores_oa(fib_stiff, partcoll_bools, nodalprop_bools, orient_bools, feas_bools, constrad_read, repeat3x3_bool, sidenum, n_runs, n_pop)
    feas_struct_all = struct;
    conn_struct_all = struct;
    stiffrat_struct_all = struct;
    partcoll_struct_all = struct;
    nodalprop_struct_all = struct;
    orient_struct_all = struct;
    
    n_total_members = nchoosek(sidenum^2,2);
    
    filepath = "C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\Opt Prob 2\\";
    methods = ["Int Pen","AOS","Bias Init","ACH"];
    
    if constrad_read
        filepath_constrad = "Constant Radii\\";
    else
        filepath_constrad = "Variable Radii\\";
    end
    
    if repeat3x3_bool
        filepath_repeatcase = "Repeat 3x3\\";
    else
        filepath_repeatcase = "Repeat 1x1\\";
    end
    
    if (partcoll_bools(2) || nodalprop_bools(2) || orient_bools(2))
        filename = "AOSMOEA_emoea_";
    else
        filename = "EpsilonMOEA_emoea_";
    end
    filepath2 = '';
    filename2 = '';
    constr_count = 0;
    for i = 1:4
        if (partcoll_bools(i) && nodalprop_bools(i) && ~orient_bools(i))
            constrained = strcat(methods(i)," - Partcoll and NodalProp\\");
            filename_snippet = strcat("pncon",num2str(i-1),"_");
        elseif (partcoll_bools(i) && ~nodalprop_bools(i) && ~orient_bools(i))
            constrained = strcat(methods(i)," - Partcoll\\");
            filename_snippet = strcat("pcon",num2str(i-1),"_");
        elseif (~partcoll_bools(i) && nodalprop_bools(i) && ~orient_bools(i))
            constrained = strcat(methods(i)," - NodalProp\\");
            filename_snippet = strcat("ncon",num2str(i-1),"_");
        elseif (~partcoll_bools(i) && ~nodalprop_bools(i) && orient_bools(i))
            constrained = strcat(methods(i)," - Orient\\");
            filename_snippet = strcat("ocon",num2str(i-1),"_");
        elseif (partcoll_bools(i) && ~nodalprop_bools(i) && orient_bools(i))
            constrained = strcat(methods(i)," - Partcoll and Orient\\");
            filename_snippet = strcat("pocon",num2str(i-1),"_");
        elseif (~partcoll_bools(i) && nodalprop_bools(i) && orient_bools(i))
            constrained = strcat(methods(i)," - NodalProp and Orient\\");
            filename_snippet = strcat("nocon",num2str(i-1),"_");
        elseif (partcoll_bools(i) && nodalprop_bools(i) && orient_bools(i))
            constrained = strcat(methods(i)," - Partcoll Nodalprop and Orient\\");
            filename_snippet = strcat("pnocon",num2str(i-1),"_");
        else
            constrained = '';
            filename_snippet = '';
            constr_count = constr_count + 1;
        end
        filepath2 = strcat(filepath2, constrained);
        filename2 = strcat(filename2, filename_snippet);
    end
    
    filepath_feas = '';
    filename_feas = '';
    if (feas_bools(1) && ~feas_bools(2))
        filepath_feas = 'Bias Init - Inters\\';
        filename_feas = 'feasbias_';
    elseif (~feas_bools(1) && feas_bools(2))
        filepath_feas = 'AOS - Inters\\';
        filename_feas = 'feasaos_';
    elseif (feas_bools(1) && feas_bools(2))
        filepath_feas = 'Bias Init AOS - Inters\\';
        filename_feas = 'feasbiasaos_';
    end
    filepath2 = strcat(filepath2,filepath_feas);
    filename2 = strcat(filename2,filename_feas);
    
    filepath_moea = '';
    if (constr_count == 4)
        filepath_moea = "Epsilon MOEA\\";
    end
    
    if fib_stiff
        filepath3 = "Fibre Stiffness\\";
        if constrad_read
            filename2 = strcat(filename2,"_prob2_fibre.csv");
        else 
            filename2 = strcat(filename2,"_fibre_varrad.csv");
        end
    else
        filepath3 = "Truss Stiffness\\";
        if constrad_read
            filename2 = strcat(filename2,"_prob2_truss.csv");
        else
            filename2 = strcat(filename2,"_truss_varrad.csv");
        end
    end
      
    %%%% read appropriate files 
    for i = 1:n_runs
        run_num = i-1;
        
        full_filepath = strcat(filepath,filepath_repeatcase,filepath_constrad,filepath2,filepath_moea,filepath3,filename,num2str(run_num),filename2);
           
        if constrad_read
            format_string = '%s';
            for j = 1:10
                format_string = strcat(format_string,'%f');
            end
            data_table = readtable(full_filepath,'Format',format_string,'HeaderLines',1);
        else
            format_string = '';
            for j = 1:(n_total_members+10)
                format_string = strcat(format_string,'%f');
            end
            data_table = readtable(full_filepath,'Format',format_string,'HeaderLines',1);
        end

        %%%% store retrieved data into different variables
        %%%% csv_data includes: [Pen. Obj. 1, Pen.Obj. 2, True Obj. 1, True Obj. 2, Feasibility Score,
        %%%% Connectivity Score, Stiffness Ratio Constraint, Partial Collapsibility Score, 
        %%%% Nodal Properties Score, Orientation Score]
        pop_size =  size(data_table,1);
        csv_data = zeros(pop_size,10);
        if constrad_read
            designs = strings(pop_size);
            designs = data_table(:,1);
        else
            designs = zeros(pop_size,n_total_members);
            designs = data_table(:,1:n_total_members);
        end
        csv_data = data_table(:,end-9:end);
        
        csv_data_array = table2array(csv_data);
        designs_array = table2array(designs);
        
        data_array_nonans_bool = any(isnan(csv_data_array),2);
        data_array_nonans = csv_data_array(~data_array_nonans_bool,:);
    
        f_penalized = data_array_nonans(:,1:2);
        feas_array = data_array_nonans(:,5);
        conn_array = data_array_nonans(:,6);
        if constrad_read
            stiffrat_array = data_array_nonans(:,7);
        else
            if repeat3x3_bool
                stiffrat_array = data_array_nonans(:,7);
            else
                stiffrat_array = 1 - data_array_nonans(:,7); % current results csv file stores (1 - stiffrat) instead of stiffrat, rectified in java code in post 
            end
        end
        partcoll_array = data_array_nonans(:,8);
        nodalprop_array = data_array_nonans(:,9);
        orient_array = data_array_nonans(:,10);
        
        pareto_bool = paretofront(f_penalized);
        designs_pareto = designs_array(pareto_bool==1);
        feas_pareto = feas_array(pareto_bool==1);
        conn_pareto = conn_array(pareto_bool==1);
        stiffrat_pareto = stiffrat_array(pareto_bool==1);
        partcoll_pareto = partcoll_array(pareto_bool==1);
        nodalprop_pareto = nodalprop_array(pareto_bool==1); 
        orient_pareto = orient_array(pareto_bool==1);
        
        current_field = strcat('run_',num2str(i));
        
        feas_struct_all.(current_field) = feas_pareto;
        conn_struct_all.(current_field) = conn_pareto;
        stiffrat_struct_all.(current_field) = stiffrat_pareto;
        partcoll_struct_all.(current_field) = partcoll_pareto;
        nodalprop_struct_all.(current_field) = nodalprop_pareto;
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
