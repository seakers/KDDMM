%% Plot constraint satisfaction plots across NFEs
clear
close all
clc

%% Parameter definitions 
sidenum = 3;
n_total_members = nchoosek(sidenum^2,2);
choice_of_model = "Truss"; % "Fibre" -> fibre model, "Truss" -> truss model, "Beam" -> beam model

%% Cases to consider for GA data (change for appropriate problem and heuristic)
truss_problem = true; % true -> truss problem, false -> artery problem
constrad_read = true;

% Case 1 - Epsilon MOEA
case1_partcoll_bools = [false, false, false, false];
case1_nodalprop_bools = [false, false, false, false];
case1_orient_bools = [false, false, false, false];
case1_inters_bools = [false, false, false, false];

% Case 2 - Interior Penalty
case2_partcoll_bools = [false, false, false, false];
case2_nodalprop_bools = [false, false, false, false];
case2_orient_bools = [true, false, false, false];
case2_inters_bools = [false, false, false, false];

% Case 3 - AOS
case3_partcoll_bools = [false, false, false, false];
case3_nodalprop_bools = [false, false, false, false];
case3_orient_bools = [false, true, false, false];
case3_inters_bools = [false, false, false, false];

% Case 4 - Biased Initialization
case4_partcoll_bools = [false, false, false, false];
case4_nodalprop_bools = [false, false, false, false];
case4_orient_bools = [false, false, true, false];
case4_inters_bools = [false, false, false, false];

% Case 5 - ACH
case5_partcoll_bools = [false, false, false, false];
case5_nodalprop_bools = [false, false, false, false];
case5_orient_bools = [false, false, false, true];
case5_inters_bools = [false, false, false, false];

%% Extract constraint values for runs of different cases
nfe = 5000; % number of function evaluations
n_runs = 30; % number of runs

nfe_res = 250;
n_datapoints = (nfe/nfe_res) + 1;
nfe_array = 0:nfe_res:nfe;

feas_allconstr_allruns = zeros(n_datapoints,n_runs,5); 
% slice 1 -> emoea, slice 2 -> intpen, slice 3 -> aos, slice 4 -> biasinit,
% slice 5 -> ach. Each slice is an (n_datapoints x n_runs) array

conn_allconstr_allruns = zeros(n_datapoints,n_runs,5);
% slice 1 -> emoea, slice 2 -> intpen, slice 3 -> aos, slice 4 -> biasinit,
% slice 5 -> ach. Each slice is an (n_datapoints x n_runs) array

if truss_problem
    stiffrat_allconstr_allruns = zeros(n_datapoints,n_runs,5);
    % slice 1 -> emoea, slice 2 -> intpen, slice 3 -> aos, slice 4 -> biasinit,
    % slice 5 -> ach. Each slice is an (n_datapoints x n_runs) array
end

for i = 1:n_runs
    [nfe_array_emoea_run, constr_array_emoea_run, ~] = constr_from_csv_data(truss_problem, choice_of_model, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_inters_bools, constrad_read, n_total_members, i-1);
    [nfe_array_intpen_run, constr_array_intpen_run, ~] = constr_from_csv_data(truss_problem, choice_of_model, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_inters_bools, constrad_read, n_total_members, i-1);
    [nfe_array_aos_run, constr_array_aos_run, ~] = constr_from_csv_data(truss_problem, choice_of_model, case3_partcoll_bools, case3_nodalprop_bools, case3_orient_bools, case3_inters_bools, constrad_read, n_total_members, i-1);
    [nfe_array_biasinit_run, constr_array_biasinit_run, ~] = constr_from_csv_data(truss_problem, choice_of_model, case4_partcoll_bools, case4_nodalprop_bools, case4_orient_bools, case4_inters_bools, constrad_read, n_total_members, i-1);
    [nfe_array_ach_run, constr_array_ach_run, ~] = constr_from_csv_data(truss_problem, choice_of_model, case5_partcoll_bools, case5_nodalprop_bools, case5_orient_bools, case5_inters_bools, constrad_read, n_total_members, i-1);
    
    [nfe_sorted_emoea_run, sort_indices_emoea_run] = sort(nfe_array_emoea_run);
    [nfe_sorted_intpen_run, sort_indices_intpen_run] = sort(nfe_array_intpen_run);
    [nfe_sorted_aos_run, sort_indices_aos_run] = sort(nfe_array_aos_run);
    [nfe_sorted_biasinit_run, sort_indices_biasinit_run] = sort(nfe_array_biasinit_run);
    [nfe_sorted_ach_run, sort_indices_ach_run] = sort(nfe_array_ach_run);
    
    constr_sorted_emoea_run = constr_array_emoea_run(sort_indices_emoea_run,:);
    constr_sorted_intpen_run = constr_array_intpen_run(sort_indices_intpen_run,:);
    constr_sorted_aos_run = constr_array_aos_run(sort_indices_aos_run,:);
    constr_sorted_biasinit_run = constr_array_biasinit_run(sort_indices_biasinit_run,:);
    constr_sorted_ach_run = constr_array_ach_run(sort_indices_ach_run,:);
        
    for j = 1:n_datapoints
        idx_nfe_datapoint_emoea = find(nfe_sorted_emoea_run <= nfe_array(j), 1, 'last');
        idx_nfe_datapoint_intpen = find(nfe_sorted_intpen_run <= nfe_array(j), 1, 'last');
        idx_nfe_datapoint_aos = find(nfe_sorted_aos_run <= nfe_array(j), 1, 'last');
        idx_nfe_datapoint_biasinit = find(nfe_sorted_biasinit_run <= nfe_array(j), 1, 'last');
        idx_nfe_datapoint_ach = find(nfe_sorted_ach_run <= nfe_array(j), 1, 'last');
        
        feas_allconstr_allruns(j,i,1) = constr_sorted_emoea_run(idx_nfe_datapoint_emoea,1);
        feas_allconstr_allruns(j,i,2) = constr_sorted_intpen_run(idx_nfe_datapoint_intpen,1);
        feas_allconstr_allruns(j,i,3) = constr_sorted_aos_run(idx_nfe_datapoint_aos,1);
        feas_allconstr_allruns(j,i,4) = constr_sorted_biasinit_run(idx_nfe_datapoint_biasinit,1);
        feas_allconstr_allruns(j,i,5) = constr_sorted_ach_run(idx_nfe_datapoint_ach,1);
        
        conn_allconstr_allruns(j,i,1) = constr_sorted_emoea_run(idx_nfe_datapoint_emoea,2);
        conn_allconstr_allruns(j,i,2) = constr_sorted_intpen_run(idx_nfe_datapoint_intpen,2);
        conn_allconstr_allruns(j,i,3) = constr_sorted_aos_run(idx_nfe_datapoint_aos,2);
        conn_allconstr_allruns(j,i,4) = constr_sorted_biasinit_run(idx_nfe_datapoint_biasinit,2);
        conn_allconstr_allruns(j,i,5) = constr_sorted_ach_run(idx_nfe_datapoint_ach,2);
        
        if truss_problem
            stiffrat_allconstr_allruns(j,i,1) = constr_sorted_emoea_run(idx_nfe_datapoint_emoea,3);
            stiffrat_allconstr_allruns(j,i,2) = constr_sorted_intpen_run(idx_nfe_datapoint_intpen,3);
            stiffrat_allconstr_allruns(j,i,3) = constr_sorted_aos_run(idx_nfe_datapoint_aos,3);
            stiffrat_allconstr_allruns(j,i,4) = constr_sorted_biasinit_run(idx_nfe_datapoint_biasinit,3);
            stiffrat_allconstr_allruns(j,i,5) = constr_sorted_ach_run(idx_nfe_datapoint_ach,3);
        end        
    end
    
end

%% Plot the constraint plots
label_strs = ["Eps. MOEA","Int. Pen.","AOS","Bias. Init.","ACH"];
linecolors = ['k','m','b','r','g'];
linestyles = ['-','-','-','-','-'];

% Feasibility plot
constr_str = "Feasibility";
plot_constraint_plot(feas_allconstr_allruns,nfe_array,label_strs,linecolors,linestyles,constr_str)

% Connectivity plot
constr_str = "Connectivity";
plot_constraint_plot(conn_allconstr_allruns,nfe_array,label_strs,linecolors,linestyles,constr_str)

% Stiffness Ratio plot (if truss problem)
if truss_problem
    constr_str = "Stiffness Ratio";
    plot_constraint_plot(stiffrat_allconstr_allruns,nfe_array,label_strs,linecolors,linestyles,constr_str)
end


%% Functions

function [nfe_array, constr_data_array, design_array] = constr_from_csv_data(prob_truss, model_choice, partcoll_bools, nodalprop_bools, orient_bools, inters_bools, constrad_read, n_total_members, run_num)
    filepath = "C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\";
    methods = ["Int Pen", "AOS", "Bias Init", "ACH"];
    heurs_list = ["PartColl", "NodalProp", "Orient", "Inters"];
    heurs_abbrvs_list = ["p","n","o","i"];
    heur_bools = [partcoll_bools; nodalprop_bools; orient_bools; inters_bools];
    
    if constrad_read
        filepath_constrad = "Constant Radii\\";
    else
        filepath_constrad = "Variable Radii\\";
    end
    
    if prob_truss
        filepath_prob = "Truss Problem\\";
    else
        filepath_prob = "Artery Problem\\";
    end
    
    if (any(heur_bools(:,2)))
        filename = "AOSMOEA_emoea_";
    else
        filename = "EpsilonMOEA_emoea_";
    end
    filepath2 = '';
    filename2 = '';
    
    constr_count = 0;
    for i = 1:size(heur_bools,2)
        constraints = strcat(methods(i)," - ");
        constraints_abbr = "";
        heur_count = 0;
        for j = 1:size(heur_bools,1)
            if heur_bools(j,i)
                constraints = strcat(constraints, heurs_list(j), "\\");
                constraints_abbr = strcat(constraints_abbr, heurs_abbrvs_list(j));
            else
                heur_count = heur_count + 1;
            end
        end
        if heur_count < size(heur_bools,2)
            filepath2 = strcat(filepath2, constraints);
            filename2 = strcat(filename2, constraints_abbr, "con", num2str(i-1), "_");
        else
            constr_count = constr_count + 1;
        end
    end
        
    filepath_moea = '';
    if (constr_count == size(heur_bools,2))
        filepath_moea = "Epsilon MOEA\\";
    end
    
    switch model_choice
        case "Fibre"
            filepath3 = "Fibre Model\\";
            if truss_problem
                if constrad_read
                    filename2 = strcat(filename2,"_prob2_fibre_fullpop.csv");
                else 
                    filename2 = strcat(filename2,"_fibre_varrad_fullpop.csv");
                end
            else
                disp("Fiber stiffness model not suitable for artery problem")
                exit
            end
        case "Truss"
            filepath3 = "Truss Model\\";
            if prob_truss
                if constrad_read
                    filename2 = strcat(filename2,"_prob2_truss_fullpop.csv");
                else
                    filename2 = strcat(filename2,"_truss_varrad_fullpop.csv");
                end
            else
                filename2 = strcat(filename2,"_artery_truss_fullpop.csv");
            end
        case "Beam"
            filepath3 = "Beam Model\\";
            if prob_truss
                if constrad_read
                    filename2 = strcat(filename2,"_prob2_beam_fullpop.csv");
                else
                    filename2 = strcat(filename2,"_beam_varrad_fullpop.csv");
                end
            else
                filename2 = strcat(filename2,"_artery_beam_fullpop.csv");
            end
    end
    
    %%%% read appropriate file 
    full_filepath = strcat(filepath,filepath_prob,filepath_constrad,filepath3,filepath2,filepath_moea,filename,num2str(run_num),filename2);
    
    if prob_truss
        n_data = 12;
    else
        n_data = 11;
    end
        
    if constrad_read
        format_string = '%s';
        for j = 1:n_data
            format_string = strcat(format_string,'%f');
        end
        data_table = readtable(full_filepath,'Format',format_string,'HeaderLines',1);
    else
        format_string = '';
        for j = 1:(n_total_members+n_data)
            format_string = strcat(format_string,'%f');
        end
        data_table = readtable(full_filepath,'Format',format_string,'HeaderLines',1);
    end
    
    %%%% store retrieved data into different variables
    %%%% for the truss problem:
    %%%% csv_data includes: [NFE, Pen. Obj. 1, Pen.Obj. 2, True Obj. 1, True Obj. 2, Feasibility Score,
    %%%% Connectivity Score, Stiffness Ratio Constraint, Partial Collapsibility Score, 
    %%%% Nodal Properties Score, Orientation Score, Intersection Score]
    %%%% for the artery problem:
    %%%% csv_data includes: [NFE, Pen. Obj. 1, Pen.Obj. 2, True Obj. 1, True Obj. 2, Feasibility Score,
    %%%% Connectivity Score, Partial Collapsibility Score, Nodal Properties Score, Orientation Score, Intersection Score]
    
    pop_size =  size(data_table,1);
    csv_data = zeros(pop_size,n_data);
    
    if constrad_read
        designs = strings(pop_size);
        designs = data_table(:,1);
    else
        designs = zeros(pop_size,n_total_members);
        designs = data_table(:,1:n_total_members);
    end
    csv_data = data_table(:,2:end);
    
    data_array = table2array(csv_data);
    design_array = table2array(designs);
    
    data_array_nonans = data_array(~any(isnan(data_array),2),:);
    
    if prob_truss
        constr_data_array = data_array_nonans(:,6:9); % [feasibility, connectivity, stiffness ratio]
    else
        constr_data_array = data_array_nonans(:,6:8); % [feasibility, connectivity]
    end
    
    nfe_array = data_array_nonans(:,1);
end

function [] = plot_constraint_plot(constr_arrays, nfe_vals_array, label_vals, linecolors, linestyles, constr_name)
    % constr_arrays = [constr_emoea (n_datapoints x n_runs) | constr_intpen(n_datapoints x n_runs) |
    %                   constr_aos (n_datapoints x n_runs) | constr_biasinit (n_datapoints x n_runs) | 
    %                   constr_ach (n_datapoints x n_runs)]  so it is a (n_datapoints x n_runs x 5) array
    % label_vals, linecolors and linestyles = 1 x 5 string arrays

    figure
    for i = 1:size(constr_arrays,3)
        boxplot(transpose(constr_arrays(:,:,i)),'positions',nfe_vals_array)
        hold on
        med_constr_allruns = median(transpose(constr_arrays(:,:,i)));
        plot(nfe_vals_array, med_constr_allruns, 'Color', linecolors(i), 'LineStyle', linestyles(i))
        hold on
    end
    hold off
    legend(label_vals,'Location','Best')
    xlabel('NFE')
    ylabel(strcat(constr_name,' value'))
    title(strcat(constr_name,' satisfaction for different leveraging methods'))   
end
