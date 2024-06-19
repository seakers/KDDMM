%% Get bounds of objectives for the runs of a particular case for the materials problems
clear
close all
clc

%% Read data for all runs of given case
prob_truss = false; % if true -> truss problem, if false -> artery problem
constant_radii = true; % if true -> members have constant radii, if false -> members have discrete radii choices
model = "Truss"; % "Fibre" - fibre model, "Truss" - truss model, "Beam" - beam model

num_runs = 30; % number of runs
sidenum = 3;

% partcoll_bools = [int_pen, AOS, bias_init, ACH];
% nodalprop_bools = [int_pen, AOS, bias_init, ACH];
% orient_bools = [int_pen, AOS, bias_init, ACH];
% inters_bools = [int_pen, AOS, bias_init, ACH];

% Case 1: Epsilon MOEA
case1_partcoll_bools = [false, false, false, false];
case1_nodalprop_bools = [false, false, false, false];
case1_orient_bools = [false, false, false, false];
case1_inters_bools = [false, false, false, false];

% Case 2: AOS - All Heuristics
case2_partcoll_bools = [false, true, false, false];
case2_nodalprop_bools = [false, true, false, false];
case2_orient_bools = [false, true, false, false];
case2_inters_bools = [false, true, false, false];

% Case 3: AOS - Promising Heuristics 
case3_partcoll_bools = [false, false, false, false];
if prob_truss
	case3_nodalprop_bools = [false, false, false, false];
else
	case3_nodalprop_bools = [false, false, false, false];
end
case3_orient_bools = [false, true, false, false];
case3_inters_bools = [false, true, false, false];

% Generate combined true and penalized pareto front arrays for each run
% case for different nfe thresholds

% NFE threshold = 500/250
if prob_truss
    nfe_thresh1 = 500;
else
    nfe_thresh1 = 250;
end

disp(strcat('Generating PFs for NFE = ',num2str(nfe_thresh1)))
% Case 1
[f_true_sat_pareto_combined_case1_thresh1, des_pareto_combined_case1_thresh1] = obtain_combined_pareto_data_case(prob_truss, model, constant_radii, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_inters_bools, nfe_thresh1, sidenum, num_runs); 

% Case 2
[f_true_sat_pareto_combined_case2_thresh1, des_pareto_combined_case2_thresh1] = obtain_combined_pareto_data_case(prob_truss, model, constant_radii, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_inters_bools, nfe_thresh1, sidenum, num_runs); 

% Case 3
[f_true_sat_pareto_combined_case3_thresh1, des_pareto_combined_case3_thresh1] = obtain_combined_pareto_data_case(prob_truss, model, constant_radii, case3_partcoll_bools, case3_nodalprop_bools, case3_orient_bools, case3_inters_bools, nfe_thresh1, sidenum, num_runs); 

truesatobj1_max_case1_thresh1 = max(f_true_sat_pareto_combined_case1_thresh1(:,1));
truesatobj1_max_case2_thresh1 = max(f_true_sat_pareto_combined_case2_thresh1(:,1));
truesatobj1_max_case3_thresh1 = max(f_true_sat_pareto_combined_case3_thresh1(:,1));
truesatobj1_max_thresh1 = max([truesatobj1_max_case1_thresh1, truesatobj1_max_case2_thresh1, truesatobj1_max_case3_thresh1]);

truesatobj2_min_case1_thresh1 = min(f_true_sat_pareto_combined_case1_thresh1(:,2));
truesatobj2_min_case2_thresh1 = min(f_true_sat_pareto_combined_case2_thresh1(:,2));
truesatobj2_min_case3_thresh1 = min(f_true_sat_pareto_combined_case3_thresh1(:,2));
truesatobj2_min_thresh1 = min([truesatobj2_min_case1_thresh1, truesatobj2_min_case2_thresh1, truesatobj2_min_case3_thresh1]);

utopia_truesat_thresh1 = [truesatobj1_max_thresh1, truesatobj2_min_thresh1];

% NFE threshold = 1000/500
if prob_truss
    nfe_thresh2 = 1000;
else
    nfe_thresh2 = 500;
end
disp(strcat('Generating PFs for NFE = ',num2str(nfe_thresh2)))
% Case 1
[f_true_sat_pareto_combined_case1_thresh2, des_pareto_combined_case1_thresh2] = obtain_combined_pareto_data_case(prob_truss, model, constant_radii, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_inters_bools, nfe_thresh2, sidenum, num_runs); 

% Case 2
[f_true_sat_pareto_combined_case2_thresh2, des_pareto_combined_case2_thresh2] = obtain_combined_pareto_data_case(prob_truss, model, constant_radii, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_inters_bools, nfe_thresh2, sidenum, num_runs); 

% Case 3
[f_true_sat_pareto_combined_case3_thresh2, des_pareto_combined_case3_thresh2] = obtain_combined_pareto_data_case(prob_truss, model, constant_radii, case3_partcoll_bools, case3_nodalprop_bools, case3_orient_bools, case3_inters_bools, nfe_thresh2, sidenum, num_runs); 

truesatobj1_max_case1_thresh2 = max(f_true_sat_pareto_combined_case1_thresh2(:,1));
truesatobj1_max_case2_thresh2 = max(f_true_sat_pareto_combined_case2_thresh2(:,1));
truesatobj1_max_case3_thresh2 = max(f_true_sat_pareto_combined_case3_thresh2(:,1));
truesatobj1_max_thresh2 = max([truesatobj1_max_case1_thresh2, truesatobj1_max_case2_thresh2, truesatobj1_max_case3_thresh2]);

truesatobj2_min_case1_thresh2 = min(f_true_sat_pareto_combined_case1_thresh2(:,2));
truesatobj2_min_case2_thresh2 = min(f_true_sat_pareto_combined_case2_thresh2(:,2));
truesatobj2_min_case3_thresh2 = min(f_true_sat_pareto_combined_case3_thresh2(:,2));
truesatobj2_min_thresh2 = min([truesatobj2_min_case1_thresh2, truesatobj2_min_case2_thresh2, truesatobj2_min_case3_thresh2]);

utopia_truesat_thresh2 = [truesatobj1_max_thresh2, truesatobj2_min_thresh2];	

% NFE threshold = 3000/1500
if prob_truss
    nfe_thresh3 = 3000;
else
    nfe_thresh3 = 1500;
end
disp(strcat('Generating PFs for NFE = ',num2str(nfe_thresh3)))
% Case 1
[f_true_sat_pareto_combined_case1_thresh3, des_pareto_combined_case1_thresh3] = obtain_combined_pareto_data_case(prob_truss, model, constant_radii, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_inters_bools, nfe_thresh3, sidenum, num_runs); 

% Case 2
[f_true_sat_pareto_combined_case2_thresh3, des_pareto_combined_case2_thresh3] = obtain_combined_pareto_data_case(prob_truss, model, constant_radii, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_inters_bools, nfe_thresh3, sidenum, num_runs); 

% Case 3
[f_true_sat_pareto_combined_case3_thresh3, des_pareto_combined_case3_thresh3] = obtain_combined_pareto_data_case(prob_truss, model, constant_radii, case3_partcoll_bools, case3_nodalprop_bools, case3_orient_bools, case3_inters_bools, nfe_thresh3, sidenum, num_runs); 

truesatobj1_max_case1_thresh3 = max(f_true_sat_pareto_combined_case1_thresh3(:,1));
truesatobj1_max_case2_thresh3 = max(f_true_sat_pareto_combined_case2_thresh3(:,1));
truesatobj1_max_case3_thresh3 = max(f_true_sat_pareto_combined_case3_thresh3(:,1));
truesatobj1_max_thresh3 = max([truesatobj1_max_case1_thresh3, truesatobj1_max_case2_thresh3, truesatobj1_max_case3_thresh3]);

truesatobj2_min_case1_thresh3 = min(f_true_sat_pareto_combined_case1_thresh3(:,2));
truesatobj2_min_case2_thresh3 = min(f_true_sat_pareto_combined_case2_thresh3(:,2));
truesatobj2_min_case3_thresh3 = min(f_true_sat_pareto_combined_case3_thresh3(:,2));
truesatobj2_min_thresh3 = min([truesatobj2_min_case1_thresh3, truesatobj2_min_case2_thresh3, truesatobj2_min_case3_thresh3]);

utopia_truesat_thresh3 = [truesatobj1_max_thresh3, truesatobj2_min_thresh3];

% NFE threshold = 6000/3000
if prob_truss
    nfe_thresh4 = 6000;
else
    nfe_thresh4 = 3000;
end
    
disp(strcat('Generating PFs for NFE = ',num2str(nfe_thresh4)))
% Case 1
[f_true_sat_pareto_combined_case1_thresh4, des_pareto_combined_case1_thresh4] = obtain_combined_pareto_data_case(prob_truss, model, constant_radii, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_inters_bools, nfe_thresh4, sidenum, num_runs); 

% Case 2
[f_true_sat_pareto_combined_case2_thresh4, des_pareto_combined_case2_thresh4] = obtain_combined_pareto_data_case(prob_truss, model, constant_radii, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_inters_bools, nfe_thresh4, sidenum, num_runs); 

% Case 3
[f_true_sat_pareto_combined_case3_thresh4, des_pareto_combined_case3_thresh4] = obtain_combined_pareto_data_case(prob_truss, model, constant_radii, case3_partcoll_bools, case3_nodalprop_bools, case3_orient_bools, case3_inters_bools, nfe_thresh4, sidenum, num_runs); 

truesatobj1_max_case1_thresh4 = max(f_true_sat_pareto_combined_case1_thresh4(:,1));
truesatobj1_max_case2_thresh4 = max(f_true_sat_pareto_combined_case2_thresh4(:,1));
truesatobj1_max_case3_thresh4 = max(f_true_sat_pareto_combined_case3_thresh4(:,1));
truesatobj1_max_thresh4 = max([truesatobj1_max_case1_thresh4, truesatobj1_max_case2_thresh4, truesatobj1_max_case3_thresh4]);

truesatobj2_min_case1_thresh4 = min(f_true_sat_pareto_combined_case1_thresh4(:,2));
truesatobj2_min_case2_thresh4 = min(f_true_sat_pareto_combined_case2_thresh4(:,2));
truesatobj2_min_case3_thresh4 = min(f_true_sat_pareto_combined_case3_thresh4(:,2));
truesatobj2_min_thresh4 = min([truesatobj2_min_case1_thresh4, truesatobj2_min_case2_thresh4, truesatobj2_min_case3_thresh4]);

utopia_truesat_thresh4 = [truesatobj1_max_thresh4, truesatobj2_min_thresh4];

% Plotting fully satisfying designs in true objectives space
figure
subplot(2,2,1)
scatter(f_true_sat_pareto_combined_case1_thresh1(:,1), f_true_sat_pareto_combined_case1_thresh1(:,2), 50, 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_sat_pareto_combined_case2_thresh1(:,1), f_true_sat_pareto_combined_case2_thresh1(:,2), 50, 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_sat_pareto_combined_case3_thresh1(:,1), f_true_sat_pareto_combined_case3_thresh1(:,2), 50, 'Marker', 's', 'MarkerEdgeColor', 'black') 
% hold on
% scatter(utopia_truesat_thresh1(1), utopia_truesat_thresh1(2), 'Marker', 'p', 'MarkerEdgeColor', 'red')
hold off
ax = gca;
ax.FontSize = 16;
if prob_truss
	xlabel('$C_{22}$','Interpreter','Latex','FontSize',16)
	ylabel('$v_f$','Interpreter','Latex','FontSize',16)
else
	xlabel('$\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
	ylabel('$deviation$','Interpreter','Latex','FontSize',16)
end
% xlim([0,E/2])
% ylim([0,1])1
%legend('Eps. MOEA','All Heurs','Prom Heurs','Location','best')
title(strcat(num2str(nfe_thresh1),' NFE'),'FontSize',16)

subplot(2,2,2)
scatter(f_true_sat_pareto_combined_case1_thresh2(:,1), f_true_sat_pareto_combined_case1_thresh2(:,2), 50, 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_sat_pareto_combined_case2_thresh2(:,1), f_true_sat_pareto_combined_case2_thresh2(:,2), 50, 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_sat_pareto_combined_case3_thresh2(:,1), f_true_sat_pareto_combined_case3_thresh2(:,2), 50, 'Marker', 's', 'MarkerEdgeColor', 'black') 
% hold on
% scatter(utopia_truesat_thresh2(1), utopia_truesat_thresh2(2), 'Marker', 'p', 'MarkerEdgeColor', 'red')
hold off
ax = gca;
ax.FontSize = 16;
if prob_truss
	xlabel('$C_{22}$','Interpreter','Latex','FontSize',16)
	ylabel('$v_f$','Interpreter','Latex','FontSize',16)
else
	xlabel('$\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
	ylabel('$deviation$','Interpreter','Latex','FontSize',16)
end
% xlim([0,E/2])
% ylim([0,1])
%legend('Eps. MOEA','All Heurs','Prom Heurs','Location','best')
title(strcat(num2str(nfe_thresh2),' NFE'),'FontSize',16)

subplot(2,2,3)
scatter(f_true_sat_pareto_combined_case1_thresh3(:,1), f_true_sat_pareto_combined_case1_thresh3(:,2), 50, 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_sat_pareto_combined_case2_thresh3(:,1), f_true_sat_pareto_combined_case2_thresh3(:,2), 50, 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_sat_pareto_combined_case3_thresh3(:,1), f_true_sat_pareto_combined_case3_thresh3(:,2), 50, 'Marker', 's', 'MarkerEdgeColor', 'black') 
% hold on
% scatter(utopia_truesat_thresh3(1), utopia_truesat_thresh3(2), 'Marker', 'p', 'MarkerEdgeColor', 'red')
hold off
ax = gca;
ax.FontSize = 16;
if prob_truss
	xlabel('$C_{22}$','Interpreter','Latex','FontSize',16)
	ylabel('$v_f$','Interpreter','Latex','FontSize',16)
else
	xlabel('$\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
	ylabel('$deviation$','Interpreter','Latex','FontSize',16)
end
% xlim([0,E/2])
% ylim([0,1])
%legend('Eps. MOEA','All Heurs','Prom Heurs','Location','best')
title(strcat(num2str(nfe_thresh3),' NFE'),'FontSize',16)

subplot(2,2,4)
scatter(f_true_sat_pareto_combined_case1_thresh4(:,1), f_true_sat_pareto_combined_case1_thresh4(:,2), 50, 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_sat_pareto_combined_case2_thresh4(:,1), f_true_sat_pareto_combined_case2_thresh4(:,2), 50, 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_sat_pareto_combined_case3_thresh4(:,1), f_true_sat_pareto_combined_case3_thresh4(:,2), 50, 'Marker', 's', 'MarkerEdgeColor', 'black') 
% hold on
% scatter(utopia_truesat_thresh4(1), utopia_truesat_thresh4(2), 'Marker', 'p', 'MarkerEdgeColor', 'red')
hold off
ax = gca;
ax.FontSize = 16;
if prob_truss
	xlabel('$C_{22}$','Interpreter','Latex','FontSize',16)
	ylabel('$v_f$','Interpreter','Latex','FontSize',16)
else
	xlabel('$\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
	ylabel('$deviation$','Interpreter','Latex','FontSize',16)
end
% xlim([0,E/2])
% ylim([0,1])
%legend('Eps. MOEA','All Heurs','Prom Heurs','Location','best')
title(strcat(num2str(nfe_thresh4),' NFE'),16)
%saveas(gcf,'pareto_truefeas_nfequadchart.png')

%% Functions

function [objs_true_sat_pareto_combined, designs_true_sat_pareto_combined] = obtain_combined_pareto_data_case(truss_prob, model_choice, read_constrad, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_threshold, sidenodenum, n_runs)
    
    n_members_total = nchoosek(sidenodenum^2,2);    
    data_array_struct = struct;
    designs_array_struct = struct;
    for i = 1:n_runs
        [data_array_case, designs_array_case] = read_csv_data_tillnfe(truss_prob, model_choice, read_constrad, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_threshold, n_members_total, i-1);
        current_field = strcat('run_',num2str(i));
        data_array_struct.(current_field) = data_array_case;
        designs_array_struct.(current_field) = designs_array_case;
    end

    [true_obj1_combined, true_obj2_combined, constr_aggr_combined, des_combined] = create_combined_arrays(data_array_struct, designs_array_struct, read_constrad, truss_prob, n_members_total, n_runs);
    f_true_combined = [-true_obj1_combined, true_obj2_combined];
    pareto_bool = compute_pareto_front_constrained(f_true_combined, constr_aggr_combined);
    
    objs_true_sat_pareto_combined = [true_obj1_combined(pareto_bool==1 & constr_aggr_combined==0), true_obj2_combined(pareto_bool==1 & constr_aggr_combined==0)];
    
    designs_true_sat_pareto_combined = des_combined(pareto_bool==1 & constr_aggr_combined==0,:);
    	
end

function [data_array_req, design_array_req] = read_csv_data_tillnfe(problem_truss, choice_of_model, constrad_read, partcoll_bools, nodalprop_bools, orient_bools, inters_bools, nfe_to_reach, n_total_members, run_num)
    filepath = "C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\"; % for lab system 
    %filepath = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\KD3M3\\result\\"; % for home system 
    methods = ["Int Pen", "AOS", "Bias Init", "ACH"];
    heurs_list = ["PartColl", "NodalProp", "Orient", "Inters"];
    heurs_abbrvs_list = ["p","n","o","i"];
    heur_bools = [partcoll_bools; nodalprop_bools; orient_bools; inters_bools];
	
    if constrad_read
        filepath_constrad = "Constant Radii\\";
    else
        filepath_constrad = "Variable Radii\\";
    end
    
    if problem_truss
        filepath_prob = "Truss Problem\\";
    else
        filepath_prob = "Artery Problem\\";
    end
    
    if (any(heur_bools(:,2)))
        filename = "AOSMOEA_emoea_";
    else
        filename = "EpsilonMOEA_emoea_";
    end
    
    filepath_cred = "";
    if (any(heur_bools(:,2)))
        filepath_cred = "set contribution dominance\\";
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
                constraints = strcat(constraints, heurs_list(j));
                constraints_abbr = strcat(constraints_abbr, heurs_abbrvs_list(j));
            else
                heur_count = heur_count + 1;
            end
        end
        if heur_count < size(heur_bools,2)
            filepath2 = strcat(filepath2, constraints, "\\");
            filename2 = strcat(filename2, constraints_abbr, "con", num2str(i-1), "_");
        else
            constr_count = constr_count + 1;
        end
    end
    
    filepath_moea = '';
    if (constr_count == size(heur_bools,2))
        filepath_moea = "Epsilon MOEA\\";
    end
    
    switch choice_of_model
        case "Fibre"
            filepath3 = "Fibre Model\\";
            if problem_truss
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
            if problem_truss
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
            if problem_truss
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
    full_filepath = strcat(filepath,filepath_prob,filepath_constrad,filepath3,filepath2,filepath_moea,filepath_cred,filename,num2str(run_num),filename2);
    
    if problem_truss
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
    %%%% Nodal Properties Score, Orientation Score]
    %%%% for the artery problem:
    %%%% csv_data includes: [NFE, Pen. Obj. 1, Pen.Obj. 2, True Obj. 1, True Obj. 2, Feasibility Score,
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
    
    full_nfe_array = data_array(:,1);
    [nfe_sorted, sort_indices] = sort(full_nfe_array);
    data_array_sorted = data_array(sort_indices,:);
    design_array_sorted = design_array(sort_indices,:);
    
    [~,closest_nfe_index] = min(abs(nfe_sorted - nfe_to_reach));
    %nfe_array = nfe_sorted(1:closest_nfe_index);
    data_array_req = data_array_sorted(1:closest_nfe_index,:);
    design_array_req = design_array_sorted(1:closest_nfe_index,:);
    
end

function [true_obj1_combined, true_obj2_combined, aggr_constr_combined, designs_combined] = create_combined_arrays(data_arrays_struct, designs_arrays_struct, read_constrad, truss_prob, n_total_members, n_runs)
    n_total = 0;
    data_array_run0 = data_arrays_struct.('run_1');
    n_constr = size(data_array_run0,2) - 4 - 4 - 1; % number of constraints changes based on problem, so subtract 4 (2 pen. and 2 true objs.), 4 (heurs) and 1 (NFE) from number of total columns
    for i = 1:n_runs
        current_field = strcat('run_',num2str(i));
        current_data_array = data_arrays_struct.(current_field);
        data_array_nonans_bool = any(isnan(current_data_array),2);
        data_array_nonans = current_data_array(~data_array_nonans_bool,:);
        n_total = n_total + size(data_array_nonans(:,1),1);
    end
    true_obj1_combined = zeros(n_total,1);
    true_obj2_combined = zeros(n_total,1);
    aggr_constr_combined = zeros(n_total,1);
    if read_constrad
        designs_combined = strings(n_total,1);
    else
        designs_combined = zeros(n_total,n_total_members);
    end
    index = 1;
    for i = 1:n_runs
        current_field = strcat('run_',num2str(i));
        current_data_array = data_arrays_struct.(current_field);
        data_array_nonans_bool = any(isnan(current_data_array),2);
        data_array_nonans = current_data_array(~data_array_nonans_bool,:);
        n_current = size(data_array_nonans(:,1),1);
        true_obj1_combined(index:index+n_current-1,1) = data_array_nonans(:,4);
        true_obj2_combined(index:index+n_current-1,1) = data_array_nonans(:,5);
        constr_combined = data_array_nonans(:,6:6+n_constr-1);
		for j = 1:size(constr_combined,1)
            feas_violation = 1 - constr_combined(j,1);
            conn_violation = 1 - constr_combined(j,2);
            if truss_prob
                stiffrat_violation = constr_combined(j,3);
                constr_violations = [feas_violation, conn_violation, stiffrat_violation];
            else
                constr_violations = [feas_violation, conn_violation];
            end
			aggr_constr_combined(index+j-1,1) = sum(abs(constr_violations)); 
		end
        if read_constrad
            current_designs_array = designs_arrays_struct.(current_field);
            design_array_nonans = current_designs_array(~data_array_nonans_bool,:);
            designs_combined(index:index+n_current-1,1) = design_array_nonans(:,1);
        else
            current_designs_array = designs_arrays_struct.(current_field);
            design_array_nonans = current_designs_array(~data_array_nonans_bool,:);
            designs_combined(index:index+n_current-1,:) = design_array_nonans;
        end
             
        index = index + n_current;
    end
end