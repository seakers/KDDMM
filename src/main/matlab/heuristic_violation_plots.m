%% Heuristic Violation Plots (truss and artery problems) - Truss and Beam Models (only GA runs)
clear 
close all
clc

%% Cases to consider for GA data
constrad_read = true;
% Case - Epsilon MOEA
truss_problem = true; % true -> truss problem, false -> artery problem
case1_partcoll_bools = [false, false, false, false];
case1_nodalprop_bools = [false, false, false, false];
case1_orient_bools = [false, false, false, false];
case1_inters_bools = [false, false, false, false];

c_ratio = 1;
if ~truss_problem
    c_ratio = 0.421;
end

choice_of_model = "Truss"; % "Fibre" -> fibre model, "Truss" -> truss model, "Beam" -> beam model
sidenum = 3;

%% Extract GA data for plotting
n_runs = 10; % number of runs to generate "n_des" architectures

min_dist_pen_pf_all = struct;
min_dist_true_pf_all = struct;

obj1_true_all = struct;
obj2_true_all = struct;
obj1_norm_all = struct;
obj2_norm_all = struct;

feas_all = struct;
conn_all = struct;
if truss_problem
    stiff_rat_all = struct;
end

coll_all = struct;
nod_all = struct;
orient_all = struct;
inters_all = struct;

[f_norm_nonans_allcases, f_true_nonans_allcases, constr_nonans_allcases, heur_nonans_allcases, ~] = obtain_combined_data_allruns(truss_problem, choice_of_model, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_inters_bools, constrad_read, sidenum, n_runs);   

for i = 1:n_runs   
    current_field = strcat('trial',num2str(i));
    f_norm_nonans_currentcase = f_norm_nonans_allcases.(current_field);
    f_true_nonans_currentcase = f_true_nonans_allcases.(current_field);
    constr_nonans_currentcase = constr_nonans_allcases.(current_field);
    heur_nonans_currentcase = heur_nonans_allcases.(current_field);

    % Penalized objectives
    obj1_norm_total = f_norm_nonans_currentcase(:,1);
    obj2_norm_total = f_norm_nonans_currentcase(:,2);

    % True objectives
    obj1_true_total = f_true_nonans_currentcase(:,1);
    obj2_true_total = f_true_nonans_currentcase(:,2);

    % Constraints
    feas_total = constr_nonans_currentcase(:,1);
    conn_total = constr_nonans_currentcase(:,2);
    if truss_problem
        stiffrat_total = constr_nonans_currentcase(:,3);
    end

    % Heuristics
    coll_total = heur_nonans_currentcase(:,1);
    nod_total = heur_nonans_currentcase(:,2);
    orient_total = heur_nonans_currentcase(:,3);
    inters_total = heur_nonans_currentcase(:,4);

    feas_all.(current_field) = feas_total;
    conn_all.(current_field) = conn_total;
    if truss_problem
        stiff_rat_all.(current_field) = stiffrat_total;
    end
    
    obj1_true_all.(current_field) = obj1_true_total;
    obj2_true_all.(current_field) = obj2_true_total;
    
    obj1_norm_all.(current_field) = obj1_norm_total;
    obj2_norm_all.(current_field) = obj2_norm_total;

    % Heuristics
    coll_all.(current_field) = coll_total;
    nod_all.(current_field) = nod_total;
    orient_all.(current_field) = orient_total;
    inters_all.(current_field) = inters_total;

end

%% Heuristic violation plots

n_des_total = 0;
for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    obj1_total = obj1_true_all.(current_field);
    n_des_total = n_des_total + size(obj1_total,1);
end

% Objectives
obj1_norm_array = zeros(n_des_total,1);
obj2_norm_array = zeros(n_des_total,1);

obj1_true_array = zeros(n_des_total,1);
obj2_true_array = zeros(n_des_total,1);

% Constraints
feas_array = zeros(n_des_total,1);
conn_array = zeros(n_des_total,1);
if truss_problem
    stiffrat_array = zeros(n_des_total,1);
end

% Heuristics
coll_array = zeros(n_des_total,1);
nod_array = zeros(n_des_total,1);
orient_array = zeros(n_des_total,1);
inters_array = zeros(n_des_total,1);

index = 1;
for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    % Objectives
    obj1_true_total = obj1_true_all.(current_field);
    obj2_true_total = obj2_true_all.(current_field);
    
    obj1_norm_total = obj1_norm_all.(current_field);
    obj2_norm_total = obj2_norm_all.(current_field);
    
    % Constraints
    feas_total = feas_all.(current_field);
    conn_total = conn_all.(current_field);
    if truss_problem
        stiffrat_total = stiff_rat_all.(current_field);
    end
    
    % Heuristics
    coll_total = coll_all.(current_field);
    nod_total = nod_all.(current_field);
    orient_total = orient_all.(current_field);
    inters_total = inters_all.(current_field);

    n_des_run = size(obj1_true_total,1);

    % Populate Objectives
    obj1_true_array(index:index+n_des_run-1,1) = obj1_true_total;
    obj2_true_array(index:index+n_des_run-1,1) = obj2_true_total;
    
    obj1_norm_array(index:index+n_des_run-1,1) = obj1_norm_total;
    obj2_norm_array(index:index+n_des_run-1,1) = obj2_norm_total;
    
    % Populate Constraints
    feas_array(index:index+n_des_run-1,1) = feas_total;
    conn_array(index:index+n_des_run-1,1) = conn_total;
    if truss_problem
        stiffrat_array(index:index+n_des_run-1,1) = stiffrat_total;
    end
    
    % Popoulate Heuristics
    coll_array(index:index+n_des_run-1,1) = coll_total;
    nod_array(index:index+n_des_run-1,1) = nod_total;
    orient_array(index:index+n_des_run-1,1) = orient_total;
    inters_array(index:index+n_des_run-1,1) = inters_total;

    index = index + n_des_run;
end

figure 
scatter(obj1_norm_array,obj2_norm_array,[],1 - coll_array,'filled')
if truss_problem
    xlabel('Normalized $C_{22}$','Interpreter','Latex','FontSize',16)
    ylabel('Normalized $v_f$','Interpreter','Latex','FontSize',16)
else 
    xlabel('Normalized $\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
    ylabel('Normalized deviation','Interpreter','Latex','FontSize',16)
end
colorbar
title('Partial Collapsibility Violation','FontSize',16)

figure 
scatter(obj1_norm_array,obj2_norm_array,[],1 - nod_array,'filled')
if truss_problem
    xlabel('Normalized $C_{22}$','Interpreter','Latex','FontSize',16)
    ylabel('Normalized $v_f$','Interpreter','Latex','FontSize',16)
else
    xlabel('Normalized $\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
    ylabel('Normalized deviation','Interpreter','Latex','FontSize',16)
end
colorbar;
title('Nodal Properties Violation','FontSize',16)

figure 
scatter(obj1_norm_array,obj2_norm_array,[],1 - orient_array,'filled')
if truss_problem
    xlabel('Normalized $C_{22}$','Interpreter','Latex','FontSize',16)
    ylabel('Normalized $v_f$','Interpreter','Latex','FontSize',16)
else
    xlabel('Normalized $\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
    ylabel('Normalized deviation','Interpreter','Latex','FontSize',16)
end
colorbar;
title('Orientation Violation','FontSize',16)

figure 
scatter(obj1_norm_array,obj2_norm_array,[],1 - inters_array,'filled')
if truss_problem
    xlabel('Normalized $C_{22}$','Interpreter','Latex','FontSize',16)
    ylabel('Normalized $v_f$','Interpreter','Latex','FontSize',16)
else
    xlabel('Normalized $\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
    ylabel('Normalized deviation','Interpreter','Latex','FontSize',16)
end
colorbar;
title('Intersection Violation','FontSize',16)

figure 
scatter(obj1_true_array,obj2_true_array,[],1 - coll_array,'filled')
if truss_problem
    xlabel('$C_{22}$','Interpreter','Latex','FontSize',16)
    ylabel('$v_f$','Interpreter','Latex','FontSize',16)
else
    xlabel('$\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
    ylabel('deviation','Interpreter','Latex','FontSize',16)
end
colorbar
title('Partial Collapsibility Violation','FontSize',16)

figure 
scatter(obj1_true_array,obj2_true_array,[],1 - nod_array,'filled')
if truss_problem
    xlabel('$C_{22}$','Interpreter','Latex','FontSize',16)
    ylabel('$v_f$','Interpreter','Latex','FontSize',16)
else
    xlabel('$\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
    ylabel('deviation','Interpreter','Latex','FontSize',16)
end
colorbar;
title('Nodal Properties Violation','FontSize',16)

figure 
scatter(obj1_true_array,obj2_true_array,[],1 - orient_array,'filled')
if truss_problem
    xlabel('$C_{22}$','Interpreter','Latex','FontSize',16)
    ylabel('$v_f$','Interpreter','Latex','FontSize',16)
else
    xlabel('$\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
    ylabel('deviation','Interpreter','Latex','FontSize',16)
end
colorbar;
title('Orientation Violation','FontSize',16)

figure 
scatter(obj1_true_array,obj2_true_array,[],1 - inters_array,'filled')
if truss_problem
    xlabel('$C_{22}$','Interpreter','Latex','FontSize',16)
    ylabel('$v_f$','Interpreter','Latex','FontSize',16)
else
    xlabel('$\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
    ylabel('deviation','Interpreter','Latex','FontSize',16)
end
colorbar;
title('Intersection Violation','FontSize',16)

%% Functions

function [objs_norm_nonans_allcases, objs_true_nonans_allcases, constraints_nonans_allcases, heuristics_nonans_allcases, designs_nonans_allcases] = obtain_combined_data_allruns(truss_prob, model_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, constrad_prob_read, sidenodenum, n_runs)
    
    n_members_total = nchoosek(sidenodenum^2,2);    
   
    objs_norm_nonans_allcases = struct;
    objs_true_nonans_allcases = struct;
    constraints_nonans_allcases = struct;
    heuristics_nonans_allcases = struct;
    designs_nonans_allcases = struct;
    
    for i = 1:n_runs
        [data_array, designs_array] = read_csv_data(truss_prob, model_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, constrad_prob_read, n_members_total, i-1);
        
        n_constr = size(data_array,2) - 4 - 4 - 1; % number of constraints changes based on problem, so subtract 4 (2 pen. and 2 true objs.), 1 for NFE and 4 (heurs) from number of total columns
        data_array_nonans_bool = any(isnan(data_array),2);
        data_array_nonans = data_array(~data_array_nonans_bool,:);
        designs_array_nonans = designs_array(~data_array_nonans_bool,:);
        
        current_field = strcat('trial',num2str(i));
        
        objs_norm_nonans_allcases.(current_field) = [data_array_nonans(:,2), data_array_nonans(:,3)];
        objs_true_nonans_allcases.(current_field) = [data_array_nonans(:,4), data_array_nonans(:,5)];
        constraints_nonans_allcases.(current_field) = data_array_nonans(:,6:6+n_constr-1);
        heuristics_nonans_allcases.(current_field) = data_array_nonans(:,6+n_constr:end);
        designs_nonans_allcases.(current_field) = designs_array_nonans;
        
    end
    
end

function [data_array, design_array] = read_csv_data(prob_truss, model_choice, partcoll_bools, nodalprop_bools, orient_bools, inters_bools, constrad_read, n_total_members, run_num)
    %filepath = "C:\\SEAK Lab\\SEAK Lab
    %Github\\KD3M3\\Truss_AOS\\result\\"; % for lab system 
    filepath = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\KD3M3\\result\\"; % for laptop
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
    
%     for i = 1:4
%         if (partcoll_bools(i) && nodalprop_bools(i) && ~orient_bools(i))
%             constrained = strcat(methods(i)," - Partcoll and NodalProp\\");
%             filename_snippet = strcat("pncon",num2str(i-1),"_");
%         elseif (partcoll_bools(i) && ~nodalprop_bools(i) && ~orient_bools(i))
%             constrained = strcat(methods(i)," - Partcoll\\");
%             filename_snippet = strcat("pcon",num2str(i-1),"_");
%         elseif (~partcoll_bools(i) && nodalprop_bools(i) && ~orient_bools(i))
%             constrained = strcat(methods(i)," - NodalProp\\");
%             filename_snippet = strcat("ncon",num2str(i-1),"_");
%         elseif (~partcoll_bools(i) && ~nodalprop_bools(i) && orient_bools(i))
%             constrained = strcat(methods(i)," - Orient\\");
%             filename_snippet = strcat("ocon",num2str(i-1),"_");
%         elseif (partcoll_bools(i) && ~nodalprop_bools(i) && orient_bools(i))
%             constrained = strcat(methods(i)," - Partcoll and Orient\\");
%             filename_snippet = strcat("pocon",num2str(i-1),"_");
%         elseif (~partcoll_bools(i) && nodalprop_bools(i) && orient_bools(i))
%             constrained = strcat(methods(i)," - NodalProp and Orient\\");
%             filename_snippet = strcat("nocon",num2str(i-1),"_");
%         elseif (partcoll_bools(i) && nodalprop_bools(i) && orient_bools(i))
%             constrained = strcat(methods(i)," - Partcoll Nodalprop and Orient\\");
%             filename_snippet = strcat("pnocon",num2str(i-1),"_");
%         else
%             constrained = '';
%             filename_snippet = '';
%             constr_count = constr_count + 1;
%         end
%         filepath2 = strcat(filepath2, constrained);
%         filename2 = strcat(filename2, filename_snippet);
%     end
    
    constr_count = 0;
    for i = 1:size(heur_bools,2)
        constraints = strcat(methods(i)," - ");
        constraints_abbr = "";
        heur_count = 0;
        for j = 1:size(heur_bools,1)
            if heur_bools(j,i)
                constraints = strcat(constraints, heurs_list(j), " ");
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
        
%     filepath_feas = '';
%     filename_feas = '';
%     if (inters_bools(1) && ~inters_bools(2))
%         filepath_feas = 'Bias Init - Inters\\';
%         filename_feas = 'feasbias_';
%     elseif (~inters_bools(1) && inters_bools(2))
%         filepath_feas = 'AOS - Inters\\';
%         filename_feas = 'feasaos_';
%     elseif (inters_bools(1) && inters_bools(2))
%         filepath_feas = 'Bias Init AOS - Inters\\';
%         filename_feas = 'feasbiasaos_';
%     end
%     filepath2 = strcat(filepath2,filepath_feas);
%     filename2 = strcat(filename2,filename_feas);
    
    filepath_moea = '';
    if (constr_count == size(heur_bools,2))
        filepath_moea = "Epsilon MOEA - Metrics\\";
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
    %%%% csv_data includes: [Pen. Obj. 1, Pen.Obj. 2, True Obj. 1, True Obj. 2, Feasibility Score,
    %%%% Connectivity Score, Stiffness Ratio Constraint, Partial Collapsibility Score, 
    %%%% Nodal Properties Score, Orientation Score]
    %%%% for the artery problem:
    %%%% csv_data includes: [Pen. Obj. 1, Pen.Obj. 2, True Obj. 1, True Obj. 2, Feasibility Score,
    %%%% Connectivity Score, Partial Collapsibility Score, Nodal Properties Score, Orientation Score]
    
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
end
