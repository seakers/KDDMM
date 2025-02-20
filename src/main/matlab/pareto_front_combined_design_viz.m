%% Plot combined pareto front with design visualization
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
c_ratio = 1;
sidenum = 3;
NC = generateNC(sel,sidenum);
CA_all = get_CA_all(sidenum);
nucFac = 1;
acceptable_rad_factor = 0.5; % parameter set in java project (smallestAcceptableRadiusFactor), change accordingly

%% Plotting 
% Case booleans
constrad_read_case = false;
repeat3x3_case = true;
case_partcoll_bools = [false, true, false, false];
case_nodalprop_bools = [false, false, false, true];
case_orient_bools = [false, true, true, false];

% Generate combined true and penalized pareto front arrays for case
[f_pen_pareto_combined, f_true_pareto_combined, constr_pareto_combined, heur_pareto_combined, des_pareto_combined] = obtain_combined_pareto_data_case(fibre_model_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, constrad_read_case, repeat3x3_case, sidenum, pop_size, num_runs); 

% Plotting penalized pareto fronts
figure
scatter(-f_pen_pareto_combined(:,1), f_pen_pareto_combined(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black', 'ButtonDownFcn', {@plot_penobjs_combined_pareto_callback,NC,CA_all,des_pareto_combined,constr_pareto_combined,heur_pareto_combined,f_pen_pareto_combined,r,acceptable_rad_factor,constrad_read_case}) 
xlabel('Penalized $C_{22}$','Interpreter','Latex')
ylabel('Penalized $v_f$','Interpreter','Latex')

% Plotting true pareto fronts
figure
scatter(f_true_pareto_combined(:,1), f_true_pareto_combined(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black', 'ButtonDownFcn', {@plot_trueobjs_combined_pareto_callback,NC,CA_all,des_pareto_combined,constr_pareto_combined,heur_pareto_combined,f_true_pareto_combined,r,acceptable_rad_factor,constrad_read_case}) 
xlabel('$C_{22}$','Interpreter','Latex')
ylabel('$v_f$','Interpreter','Latex')

%% Functions

function [objs_pen_pareto_combined, objs_true_pareto_combined, constraints_pareto_combined, heuristics_pareto_combined, designs_pareto_combined] = obtain_combined_pareto_data_case(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, constrad_prob_read, repeat3x3_bool, sidenodenum, n_pop, n_runs)
    
    n_members_total = nchoosek(sidenodenum^2,2);    
    data_array = zeros(n_pop,10,n_runs);
    if constrad_prob_read
        designs_array = strings(n_pop,n_runs);
        for i = 1:n_runs
            [data_array(:,:,i), designs_array(:,i)] = read_csv_data(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, constrad_prob_read, repeat3x3_bool, n_members_total, i-1);
        end
    else
        designs_array = zeros(n_pop,n_members_total,n_runs);
        for i = 1:n_runs
            [data_array(:,:,i), designs_array(:,:,i)] = read_csv_data(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, constrad_prob_read, repeat3x3_bool, n_members_total, i-1);
        end
    end
    [pen_obj1_combined, pen_obj2_combined, true_obj1_combined, true_obj2_combined, constraints_combined, heuristics_combined, des_combined] = create_combined_arrays(data_array, designs_array, constrad_prob_read, n_members_total, n_runs);
    f_pen_combined = [pen_obj1_combined, pen_obj2_combined];
    pareto_bool = paretofront(f_pen_combined);
    objs_pen_pareto_combined = f_pen_combined(pareto_bool==1,:);
    objs_true_pareto_combined = [true_obj1_combined(pareto_bool==1), true_obj2_combined(pareto_bool==1)];
    constraints_pareto_combined = constraints_combined(pareto_bool==1,:);
    heuristics_pareto_combined = heuristics_combined(pareto_bool==1,:);
    designs_pareto_combined = des_combined(pareto_bool==1,:);
    
end

function [data_array, design_array] = read_csv_data(fib_stiff, partcoll_bools, nodalprop_bools, orient_bools, constrad_read, repeat3x3_boolean, n_total_members, run_num)
    filepath = "C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\OA run data - optimization Problem 2\\";
    methods = ["Int Pen","AOS","Bias Init","ACH"];
    
    if constrad_read
        filepath_constrad = "Constant Radii\\";
    else
        filepath_constrad = "Variable Radii\\";
    end
    
    if repeat3x3_boolean
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
            constrained = strcat(methods(i)," - Orientation\\");
            filename_snippet = strcat("ocon",num2str(i-1),"_");
        elseif (partcoll_bools(i) && ~nodalprop_bools(i) && orient_bools(i))
            constrained = strcat(methods(i)," - Partcoll and Orientation\\");
            filename_snippet = strcat("pocon",num2str(i-1),"_");
        elseif (~partcoll_bools(i) && nodalprop_bools(i) && orient_bools(i))
            constrained = strcat(methods(i)," - NodalProp and Orientation\\");
            filename_snippet = strcat("nocon",num2str(i-1),"_");
        elseif (partcoll_bools(i) && nodalprop_bools(i) && orient_bools(i))
            constrained = strcat(methods(i)," - Partcoll Nodalprop and Orientation\\");
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

function [pen_obj1_combined, pen_obj2_combined, true_obj1_combined, true_obj2_combined, constr_combined, heur_combined, designs_combined] = create_combined_arrays(data_array, designs_array, read_constrad, n_total_members, n_runs)
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
    constr_combined = zeros(n_total,3);
    heur_combined = zeros(n_total,3);
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
        constr_combined(index:index+n_current-1,:) = data_array_nonans(:,5:7);
        heur_combined(index:index+n_current-1,:) = data_array_nonans(:,8:10);
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
