%% Plot combined pareto front comparison
clear all
close all
clc

%% CSV Read parameters
fibre_model_used = false;
sidenum = 3;

num_runs = 30; % change based on run 
pop_size = 100; % change based on run

% partcoll_bools = [int_pen, AOS, biased_init, ACH] boolean array
% nodalprop_bools = [int_pen, AOS, biased_init, ACH] boolean array
% orient_bools = [int_pen, AOS, biased_init, ACH] boolean array

%% Parameter definitions (printable designs)
E = 1.8162e6; % Young's Modulus for polymeric material (example: 1.8162 MPa for SIL material)
sel = 10e-3; % Unit square side length (NOT individual truss length) (in m)
r = 250e-6; % Radius for cross-sectional area of (assumed circular) truss members (in m)
A = pi*(r^2); % Cross-sectional area of truss member
%NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1]; 
%CA_all = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 2,3; 2,4; 2,5; 2,6; 2,7; 2,8; 2,9; 3,4; 3,5; 3,6; 3,7; 3,8; 3,9; 4,5; 4,6; 4,7; 4,8; 4,9; 5,6; 5,7; 5,8; 5,9; 6,7; 6,8; 6,9; 7,8; 7,9; 8,9];
c_ratio = 1;
sidenum = 3;
NC = generateNC(sel,sidenum);
CA_all = get_CA_all(sidenum);
nucFac = 3;

%% Comparing Combined Pareto Fronts (eps-MOEA vs AOS - Orientation)(constant radii problem)
% Case 1 - Epsilon MOEA
constrad_read_case1 = true;
repeat3x3_case1 = true;
case1_partcoll_bools = [false, false, false, false];
case1_nodalprop_bools = [false, false, false, false];
case1_orient_bools = [false, false, false, false];

% Case 2 - AOS - PartColl & Orientation, Bias Init - Orientation, ACH -
% NodalProp
constrad_read_case2 = true;
repeat3x3_case2 = true;
case2_partcoll_bools = [false, true, false, false];
case2_nodalprop_bools = [false, false, false, true];
case2_orient_bools = [false, true, true, false];

% Generate combined true and penalized pareto front arrays for each run case
% Case 1
[f_pen_pareto_combined_case1, f_true_pareto_combined_case1, des_pareto_combined_case1] = obtain_combined_pareto_data_case(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, constrad_read_case1, repeat3x3_case1, sidenum, pop_size, num_runs); 

% Case 2
[f_pen_pareto_combined_case2, f_true_pareto_combined_case2, des_pareto_combined_case2] = obtain_combined_pareto_data_case(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, constrad_read_case2, repeat3x3_case2, sidenum, pop_size, num_runs); 

% Plotting penalized pareto fronts
figure
scatter(-f_pen_pareto_combined_case1(:,1), f_pen_pareto_combined_case1(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(-f_pen_pareto_combined_case2(:,1), f_pen_pareto_combined_case2(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold off
xlabel('Penalized $C_{22}$','Interpreter','Latex')
ylabel('Penalized $v_f$','Interpreter','Latex')
legend('Without Heuristics','With Heuristics','Location','best')

% Plotting true pareto fronts
figure
scatter(f_true_pareto_combined_case1(:,1), f_true_pareto_combined_case1(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_pareto_combined_case2(:,1), f_true_pareto_combined_case2(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black')
hold off
xlabel('$C_{22}$','Interpreter','Latex')
ylabel('$v_f$','Interpreter','Latex')
legend('Without Heuristics','With Heuristics','Location','best')

%% Comparing Combined Pareto Fronts (eps-MOEA vs AOS - Orientation)(variable radii problem)
% Case 1 - Epsilon MOEA
constrad_read_case1 = false;
repeat3x3_case1 = true;
case1_partcoll_bools = [false, false, false, false];
case1_nodalprop_bools = [false, false, false, false];
case1_orient_bools = [false, false, false, false];

% Case 2 - AOS - PartColl & Orientation, Bias Init - Orientation, ACH -
% NodalProp
constrad_read_case2 = false;
repeat3x3_case2 = true;
case2_partcoll_bools = [false, true, false, false];
case2_nodalprop_bools = [false, false, false, true];
case2_orient_bools = [false, true, true, false];

% Generate combined true and penalized pareto front arrays for each run case
% Case 1
[f_pen_pareto_combined_case1, f_true_pareto_combined_case1, des_pareto_combined_case1] = obtain_combined_pareto_data_case(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, constrad_read_case1, repeat3x3_case1, sidenum, pop_size, num_runs); 

% Case 2
[f_pen_pareto_combined_case2, f_true_pareto_combined_case2, des_pareto_combined_case2] = obtain_combined_pareto_data_case(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, constrad_read_case2, repeat3x3_case2, sidenum, pop_size, num_runs); 

% Plotting penalized pareto fronts
figure
scatter(-f_pen_pareto_combined_case1(:,1), f_pen_pareto_combined_case1(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(-f_pen_pareto_combined_case2(:,1), f_pen_pareto_combined_case2(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold off
xlabel('Penalized $C_{22}$','Interpreter','Latex')
ylabel('Penalized $v_f$','Interpreter','Latex')
legend('Without Heuristics','With Heuristics','Location','best')

% Plotting true pareto fronts
figure
scatter(f_true_pareto_combined_case1(:,1), f_true_pareto_combined_case1(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_pareto_combined_case2(:,1), f_true_pareto_combined_case2(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black')
hold off
xlabel('$C_{22}$','Interpreter','Latex')
ylabel('$v_f$','Interpreter','Latex')
legend('Without Heuristics','With Heuristics','Location','best')

%% Functions 
function [objs_pen_pareto_combined, objs_true_pareto_combined, designs_pareto_combined] = obtain_combined_pareto_data_case(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, constrad_prob_read, repeat_case3x3, sidenodenum, n_pop, n_runs)
    
    n_members_total = nchoosek(sidenodenum^2,2);    
    data_array = zeros(n_pop,10,n_runs);
    if constrad_prob_read
        designs_array = strings(n_pop,n_runs);
        for i = 1:n_runs
            [data_array(:,:,i), designs_array(:,i)] = read_csv_data(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, constrad_prob_read, repeat_case3x3, n_members_total, i-1);
        end
    else
        designs_array = zeros(n_pop,n_members_total,n_runs);
        for i = 1:n_runs
            [data_array(:,:,i), designs_array(:,:,i)] = read_csv_data(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, constrad_prob_read, repeat_case3x3, n_members_total, i-1);
        end
    end
    [pen_obj1_combined, pen_obj2_combined, true_obj1_combined, true_obj2_combined, des_combined] = create_combined_arrays(data_array, designs_array, constrad_prob_read, n_members_total, n_runs);
    f_pen_combined = [pen_obj1_combined, pen_obj2_combined];
    pareto_bool = paretofront(f_pen_combined);
    objs_pen_pareto_combined = f_pen_combined(pareto_bool==1,:);
    objs_true_pareto_combined = [true_obj1_combined(pareto_bool==1), true_obj2_combined(pareto_bool==1)];
    designs_pareto_combined = des_combined(pareto_bool==1,:);
    
end

function [data_array, design_array] = read_csv_data(fib_stiff, partcoll_bools, nodalprop_bools, orient_bools, constrad_read, repeat3x3_bool, n_total_members, run_num)
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
    
    %%%% read appropriate file 
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
    
    data_array = table2array(csv_data);
    design_array = table2array(designs);
end

function [pen_obj1_combined, pen_obj2_combined, true_obj1_combined, true_obj2_combined, designs_combined] = create_combined_arrays(data_array, designs_array, read_constrad, n_total_members, n_runs)
    n_total = 0;
    for i = 1:n_runs
        current_data_array = data_array(:,:,i);
        data_array_nonans_bool = any(isnan(current_data_array),2);
        data_array_nonans = current_data_array(~data_array_nonans_bool,:);
        n_total = n_total + size(data_array_nonans(:,1),1);
    end
    pen_obj1_combined = zeros(n_total,1);
    pen_obj2_combined = zeros(n_total,1);
    true_obj1_combined = zeros(n_total,1);
    true_obj2_combined = zeros(n_total,1);
    if read_constrad
        designs_combined = strings(n_total,1);
    else
        designs_combined = zeros(n_total,n_total_members);
    end
    index = 1;
    for i = 1:n_runs
        current_data_array = data_array(:,:,i);
        data_array_nonans_bool = any(isnan(current_data_array),2);
        data_array_nonans = current_data_array(~data_array_nonans_bool,:);
        n_current = size(data_array_nonans(:,1),1);
        pen_obj1_combined(index:index+n_current-1,1) = data_array_nonans(:,1);
        pen_obj2_combined(index:index+n_current-1,1) = data_array_nonans(:,2);
        true_obj1_combined(index:index+n_current-1,1) = data_array_nonans(:,3);
        true_obj2_combined(index:index+n_current-1,1) = data_array_nonans(:,4);
        if read_constrad
            current_designs_array = designs_array(:,i);
            design_array_nonans = current_designs_array(~data_array_nonans_bool,:);
            designs_combined(index:index+n_current-1,1) = design_array_nonans(:,1);
        else
            current_designs_array = designs_array(:,:,i);
            design_array_nonans = current_designs_array(~data_array_nonans_bool,:);
            designs_combined(index:index+n_current-1,:) = design_array_nonans;
        end
             
        index = index + n_current;
    end
end

function NC = generateNC(sel,sidenum)
    notchvec = linspace(0,1,sidenum);
    NC = [];
    for i = 1:1:sidenum
        for j = 1:1:sidenum
            NC = [NC;notchvec(i),notchvec(j)];
        end
    end
    NC = sel.*NC;
end

function CA_all = get_CA_all(sidenum)
    n_member = 0;
    n_total_members = nchoosek(sidenum^2,2);
    CA_all = zeros(n_total_members,2);
    for i = 1:(sidenum^2)
        for j = i+1:(sidenum^2)
            CA_all(n_member+1,:) = [i,j];
            n_member = n_member + 1;
        end
    end
end

function [norm_objs_pen, norm_objs_true] = obtain_normalization_objs(pen_objs_combined_allcases, true_objs_combined_allcases)
    n_cases = numel(fieldnames(pen_objs_combined_allcases));
    norm_pen_obj1 = zeros(n_cases,2);
    norm_pen_obj2 = zeros(n_cases,2);
    for i = 1:n_cases
        fieldname = strcat('case',num2str(i));
        pen_objs_case = pen_objs_combined_allcases.(fieldname);
        norm_pen_obj1(i,:) = [min(pen_objs_case(:,1)),max(pen_objs_case(:,1))];
        norm_pen_obj2(i,:) = [min(pen_objs_case(:,2)),max(pen_objs_case(:,2))];
    end
    norm_objs_pen = [norm_pen_obj1, norm_pen_obj2];
    norm_true_obj1 = zeros(n_cases,2);
    norm_true_obj2 = zeros(n_cases,2);
    for i = 1:n_cases
        fieldname = strcat('case',num2str(i));
        true_objs_case = true_objs_combined_allcases.(fieldname);
        norm_true_obj1(i,:) = [min(true_objs_case(:,1)),max(true_objs_case(:,1))];
        norm_true_obj2(i,:) = [min(true_objs_case(:,2)),max(true_objs_case(:,2))];
    end
    norm_objs_true = [norm_true_obj1, norm_true_obj2];
end

function [pen_objs_case_normalized, pen_objs_overall_normalized, true_objs_case_normalized, true_objs_overall_normalized] = normalize_objs(pen_objs_allcases, true_objs_allcases, norm_pen_objs, norm_true_objs)
    % 1st obj to maximize, 2nd obj to minimize; pen_obj1 is -ve, the rest
    % are positive
    n_cases = numel(fieldnames((pen_objs_allcases)));
    norm_pen_obj1_overall_min = min(norm_pen_objs(:,1));
    norm_pen_obj1_overall_max = max(norm_pen_objs(:,2));
    norm_pen_obj2_overall_min = min(norm_pen_objs(:,3));
    norm_pen_obj2_overall_max = max(norm_pen_objs(:,4));
    norm_true_obj1_overall_min = min(norm_true_objs(:,1));
    norm_true_obj1_overall_max = max(norm_true_objs(:,2));
    norm_true_obj2_overall_min = min(norm_true_objs(:,3));
    norm_true_obj2_overall_max = max(norm_true_objs(:,4));
    
    pen_objs_case_normalized = struct;
    pen_objs_overall_normalized = struct;
    true_objs_case_normalized = struct;
    true_objs_overall_normalized = struct;
    for i = 1:n_cases
        fieldname = strcat('case',num2str(i));
        pen_objs_case = pen_objs_allcases.(fieldname);
        true_objs_case = true_objs_allcases.(fieldname);
        
        case_normalized_pen_objs_currentcase = zeros(size(pen_objs_case,1),2);
        overall_normalized_pen_objs_currentcase = zeros(size(pen_objs_case,1),2);
        case_normalized_true_objs_currentcase = zeros(size(true_objs_case,1),2);
        overall_normalized_true_objs_currentcase = zeros(size(true_objs_case,1),2);
        for j = 1:size(pen_objs_case,1)
            case_normalized_pen_objs_currentcase(j,:) = [(norm_pen_objs(i,2) - pen_objs_case(j,1))/(norm_pen_objs(i,2) - norm_pen_objs(i,1)),(norm_pen_objs(i,4) - pen_objs_case(j,2))/(norm_pen_objs(i,4) - norm_pen_objs(i,3))];
            overall_normalized_pen_objs_currentcase(j,:) = [(norm_pen_obj1_overall_max - pen_objs_case(j,1))/(norm_pen_obj1_overall_max - norm_pen_obj1_overall_min),(norm_pen_obj2_overall_max - pen_objs_case(j,2))/(norm_pen_obj2_overall_max - norm_pen_obj2_overall_min)];
            case_normalized_true_objs_currentcase(j,:) = [(true_objs_case(j,1) - norm_true_objs(i,1))/(norm_true_objs(i,2) - norm_true_objs(i,1)),(norm_true_objs(i,4) - true_objs_case(j,2))/(norm_true_objs(i,4) - norm_true_objs(i,3))];
            overall_normalized_true_objs_currentcase(j,:) = [(norm_true_obj1_overall_max - true_objs_case(j,1))/(norm_true_obj1_overall_max - norm_true_obj1_overall_min),(norm_true_obj2_overall_max - true_objs_case(j,2))/(norm_true_obj2_overall_max - norm_true_obj2_overall_min)];
        end
        field_name = strcat('case',num2str(i));
        pen_objs_case_normalized.(field_name) = case_normalized_pen_objs_currentcase;
        pen_objs_overall_normalized.(field_name) = overall_normalized_pen_objs_currentcase;
        true_objs_case_normalized.(field_name) = case_normalized_true_objs_currentcase;
        true_objs_overall_normalized.(field_name) = overall_normalized_true_objs_currentcase;
    end
end

function [utopia_dist] = compute_avg_utopia_distance(norm_objs)
    n_cases = numel(fieldnames(norm_objs));
    utopia_dist = zeros(n_cases,1);
    for i = 1:n_cases
        fieldname = strcat('case',num2str(i));
        norm_objs_case = norm_objs.(fieldname);
        utopia_dist_case = 0;
        for j = 1:size(norm_objs_case,1)
            current_norm_objs = norm_objs_case(j,:);
            current_utopia_dist = sqrt((1 - current_norm_objs(1))^2 + (1 - current_norm_objs(2))^2);
            utopia_dist_case = utopia_dist_case + current_utopia_dist;
        end
        utopia_dist(i) = utopia_dist_case/size(norm_objs_case,1);
    end
end

            
        
        
        