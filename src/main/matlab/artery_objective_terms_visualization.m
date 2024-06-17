%% Visualize the spread of designs for different objective terms at different NFE values for the artery problem
clear
close all
clc

%% Read data for all runs of Epsilon MOEA case
prob_truss = false; % if true -> truss problem, if false -> artery problem
constant_radii = true; % if true -> members have constant radii, if false -> members have discrete radii choices
model = "Truss"; % "Fibre" - fibre model, "Truss" - truss model, "Beam" - beam model

num_runs = 30; % number of runs
show_only_pareto = true; % if true -> show only designs in the combined pareto front, if false -> show all designs

% partcoll_bools = [int_pen, AOS, bias_init, ACH];
% nodalprop_bools = [int_pen, AOS, bias_init, ACH];
% orient_bools = [int_pen, AOS, bias_init, ACH];
% inters_bools = [int_pen, AOS, bias_init, ACH];

% Case: Epsilon MOEA
case_partcoll_bools = [false, false, false, false];
case_nodalprop_bools = [false, false, false, false];
case_orient_bools = [false, false, false, false];
case_inters_bools = [false, false, false, false];

%% Parameter definitions (printable designs)
E = 1.8162e6; % Young's Modulus for polymeric material (example: 1.8162 MPa for SIL material)
sel = 10e-3; % Unit square side length (NOT individual truss length) (in m)
r = 250e-6; % Radius for cross-sectional area of (assumed circular) truss members (in m)
A = pi*(r^2); % Cross-sectional area of truss member
sidenum = 3;
n_members_total = nchoosek(sidenum^2,2); 
c_target = 0.421;
nucFac = 3; % only used for fiber stiffness model

CA_all = get_CA_all(sidenum);
%NC = generateNC(sel, sidenum);

%n_members_repeated = 2*nchoosek(sidenum,2);
%n_variables = n_members_total - n_members_repeated;

%% Read designs for each run and compute first objective and artery problem constraints in the deviation objective

% NFE 1 - 200
nfe_thresh1 = 200;
disp(strcat('Computing objective terms for NFE = ',num2str(nfe_thresh1)))
[obj2_terms_thresh1, obj2_labels] = compute_obj2_terms(show_only_pareto, prob_truss, model, constant_radii, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_thresh1, CA_all, sel, sidenum, E, r, c_target, nucFac, num_runs);
plot_obj2_terms_nfe(obj2_terms_thresh1, obj2_labels, nfe_thresh1)

% NFE 2 - 400
nfe_thresh2 = 400;
disp(strcat('Computing objective terms for NFE = ',num2str(nfe_thresh2)))
[obj2_terms_thresh2, ~] = compute_obj2_terms(show_only_pareto, prob_truss, model, constant_radii, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_thresh2, CA_all, sel, sidenum, E, r, c_target, nucFac, num_runs);
plot_obj2_terms_nfe(obj2_terms_thresh2, obj2_labels, nfe_thresh2)

% NFE 3 - 600
nfe_thresh3 = 600;
disp(strcat('Computing objective terms for NFE = ',num2str(nfe_thresh3)))
[obj2_terms_thresh3, ~] = compute_obj2_terms(show_only_pareto, prob_truss, model, constant_radii, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_thresh3, CA_all, sel, sidenum, E, r, c_target, nucFac, num_runs);
plot_obj2_terms_nfe(obj2_terms_thresh3, obj2_labels, nfe_thresh3)

% NFE 4 - 800
nfe_thresh4 = 800;
disp(strcat('Computing objective terms for NFE = ',num2str(nfe_thresh4)))
[obj2_terms_thresh4, ~] = compute_obj2_terms(show_only_pareto, prob_truss, model, constant_radii, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_thresh4, CA_all, sel, sidenum, E, r, c_target, nucFac, num_runs);
plot_obj2_terms_nfe(obj2_terms_thresh4, obj2_labels, nfe_thresh4)

% NFE 5 - 1000
nfe_thresh5 = 1000;
disp(strcat('Computing objective terms for NFE = ',num2str(nfe_thresh5)))
[obj2_terms_thresh5, ~] = compute_obj2_terms(show_only_pareto, prob_truss, model, constant_radii, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_thresh5, CA_all, sel, sidenum, E, r, c_target, nucFac, num_runs);
plot_obj2_terms_nfe(obj2_terms_thresh5, obj2_labels, nfe_thresh5)

% NFE 6 - 2000
nfe_thresh6 = 2000;
disp(strcat('Computing objective terms for NFE = ',num2str(nfe_thresh6)))
[obj2_terms_thresh6, ~] = compute_obj2_terms(show_only_pareto, prob_truss, model, constant_radii, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_thresh6, CA_all, sel, sidenum, E, r, c_target, nucFac, num_runs);
plot_obj2_terms_nfe(obj2_terms_thresh6, obj2_labels, nfe_thresh6)

% NFE 7 - 3000
nfe_thresh7 = 3000;
disp(strcat('Computing objective terms for NFE = ',num2str(nfe_thresh7)))
[obj2_terms_thresh7, ~] = compute_obj2_terms(show_only_pareto, prob_truss, model, constant_radii, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_thresh7, CA_all, sel, sidenum, E, r, c_target, nucFac, num_runs);
plot_obj2_terms_nfe(obj2_terms_thresh7, obj2_labels, nfe_thresh7)

% NFE 8 - 4000
nfe_thresh8 = 4000;
disp(strcat('Computing objective terms for NFE = ',num2str(nfe_thresh8)))
[obj2_terms_thresh8, ~] = compute_obj2_terms(show_only_pareto, prob_truss, model, constant_radii, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_thresh8, CA_all, sel, sidenum, E, r, c_target, nucFac, num_runs);
plot_obj2_terms_nfe(obj2_terms_thresh8, obj2_labels, nfe_thresh8)

% NFE 9 - 5000
nfe_thresh9 = 5000;
disp(strcat('Computing objective terms for NFE = ',num2str(nfe_thresh9)))
[obj2_terms_thresh9, ~] = compute_obj2_terms(show_only_pareto, prob_truss, model, constant_radii, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_thresh9, CA_all, sel, sidenum, E, r, c_target, nucFac, num_runs);
plot_obj2_terms_nfe(obj2_terms_thresh9, obj2_labels, nfe_thresh9)

% NFE 10 - 6000
nfe_thresh10 = 6000;
disp(strcat('Computing objective terms for NFE = ',num2str(nfe_thresh10)))
[obj2_terms_thresh10, ~] = compute_obj2_terms(show_only_pareto, prob_truss, model, constant_radii, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_thresh10, CA_all, sel, sidenum, E, r, c_target, nucFac, num_runs);
plot_obj2_terms_nfe(obj2_terms_thresh10, obj2_labels, nfe_thresh10)

%% Functions

function [] = plot_obj2_terms_nfe(obj2_terms_thresh, term_labels, nfe_thresh)
    figure
    n_objs = size(obj2_terms_thresh, 2);
    for i = 1:n_objs
        subplot(n_objs/2,2,i)
        histogram(obj2_terms_thresh(:,i))
        xlabel(term_labels(i),'Interpreter','latex')    
    end
    sgtitle(strcat('NFE =',' ',num2str(nfe_thresh)))
end

function [obj2_terms_thresh, obj2_term_labels] = compute_obj2_terms(pareto_only, truss_prob, model_choice, constant_radii_used, partcoll_bools, nodalprop_bools, orient_bools, inters_bools, nfe_thresh, CA_full, sidelength, sidenodenum, E_mat, r_member, target_c, nucFactor, n_runs)
    [~, des_pareto_combined_thresh] = obtain_combined_pareto_data_case(pareto_only, truss_prob, model_choice, constant_radii_used, partcoll_bools, nodalprop_bools, orient_bools, inters_bools, nfe_thresh, sidenodenum, n_runs);

    n_des = size(des_pareto_combined_thresh,1);
    obj2_c22c11_vals_thresh = zeros(n_des, 1);
    obj2_c12c11_vals_thresh = zeros(n_des, 1);
    obj2_c21c11_vals_thresh = zeros(n_des, 1);
    obj2_c61_vals_thresh = zeros(n_des, 1);
    obj2_c62_vals_thresh = zeros(n_des, 1);
    obj2_c16_vals_thresh = zeros(n_des, 1);
    obj2_c26_vals_thresh = zeros(n_des, 1);
    obj2_c66c11_vals_thresh = zeros(n_des, 1);

    for i = 1:n_des
        x_des_current = get_binary_array_from_bitstring(des_pareto_combined_thresh{i});
        CA_des = CA_full(x_des_current == 1,:);
		rvar_des = r_member.*ones(1,size(CA_des,1));
        switch model_choice
            case "Fibre"
                if truss_prob
                    [C11,C22,~] = fiberStiffnessModel_rVar_V3(sel,rvar_des,E_mat,CA_des,nucFactor,sidenodenum);
                    C_des = zeros(6);
                    C_des(1,1) = C11;
                    C_des(2,2) = C22;
                else
                    disp("Fiber stiffness model not suitable for artery problem")
                    exit
                end
            case "Truss"
                [C_des, ~] = trussMetaCalc_NxN_1UC_rVar_AVar(sidenodenum,sidelength,rvar_des,E_mat,CA_des);
                %volfrac_des = calcVF_NxN_feasOnly(CA_des,r_member,sidelength,sidenodenum);
            case "Beam"
                C_des = Beam_2D_NxN_PBC(sidelength,sidenodenum,r_member,E_mat,CA_des);
                %volfrac_des = calcVF_NxN_feasOnly(CA_des,r_member,sidelength,sidenodenum);
        end
        obj2_c22c11_vals_thresh(i,1) = abs((C_des(2,2)/C_des(1,1)) - target_c)/6;
		obj2_c12c11_vals_thresh(i,1) = abs((C_des(1,2)/C_des(1,1)) - 0.0745);
		obj2_c21c11_vals_thresh(i,1) = abs((C_des(2,1)/C_des(1,1)) - 0.0745);
		obj2_c61_vals_thresh(i,1) = abs(C_des(3,1))/1.5e5;
		obj2_c62_vals_thresh(i,1) = abs(C_des(3,2))/1.5e5;
		obj2_c16_vals_thresh(i,1) = abs(C_des(1,3))/9e4;
		obj2_c26_vals_thresh(i,1) = abs(C_des(2,3))/9.5e4;
		obj2_c66c11_vals_thresh(i,1) = (abs((C_des(3,3)/C_des(1,1)) - 5.038) - 4.5)/0.5;

    end

    obj2_terms_thresh = [obj2_c22c11_vals_thresh, obj2_c12c11_vals_thresh, obj2_c21c11_vals_thresh, obj2_c61_vals_thresh, obj2_c62_vals_thresh, obj2_c16_vals_thresh, obj2_c26_vals_thresh, obj2_c66c11_vals_thresh];
    obj2_term_labels = [strcat("$\frac{\left|\frac{C_{22}}{C_{11}} - ",num2str(0.421),"\right|}{6}$"), "$\left|\frac{C_{12}}{C_{11}} - 0.0745\right|$", "$\left|\frac{C_{21}}{C_{11}} - 0.0745\right|$", "$\frac{\left|C_{61}\right|}{1.5\times 10^5}$", "$\frac{\left|C_{62}\right|}{1.5\times 10^5}$", "$\frac{\left|C_{16}\right|}{9\times 10^4}$", "$\frac{\left|C_{26}\right|}{9.5\times 10^4}$", "$\frac{\left|\frac{C_{66}}{C_{11}} - 5.038\right| - 4.5}{0.5}$"];

end

function [f_true_all, des_all] = obtain_combined_pareto_data_case(only_pareto, truss_prob, model_choice, read_constrad, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, nfe_threshold, sidenodenum, n_runs)
    
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

    [f_true_unique, i_fc, ~] = unique(f_true_combined, 'rows', 'stable'); 
    constr_aggr_unique = constr_aggr_combined(i_fc);
    des_unique = des_combined(i_fc);
    true_obj1_unique = true_obj1_combined(i_fc);
    true_obj2_unique = true_obj2_combined(i_fc);

    if only_pareto
        pareto_bool = compute_pareto_front_constrained(f_true_unique, constr_aggr_unique);
        f_true_all = [true_obj1_unique(pareto_bool==1), true_obj2_unique(pareto_bool==1)];
        des_all = des_unique(pareto_bool==1,:);
    else
        f_true_all = f_true_unique;
        des_all = des_unique;
    end
end

function [data_array_req, design_array_req] = read_csv_data_tillnfe(problem_truss, choice_of_model, constrad_read, partcoll_bools, nodalprop_bools, orient_bools, inters_bools, nfe_to_reach, n_total_members, run_num)
    %filepath = "C:\\SEAK Lab\\SEAK Lab
    %Github\\KD3M3\\Truss_AOS\\result\\"; % for lab system 
    filepath = "C:\\Users\\rosha\\Documents\\SEAK Lab Github\\KD3M3\\result\\"; % for home system 
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
    full_filepath = strcat(filepath,filepath_prob,filepath_constrad,filepath3,filepath2,filepath_moea,filename,num2str(run_num),filename2);
    
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

function x_des = get_binary_array_from_bitstring(des_string)
	x_des = zeros(strlength(des_string),1);
	for i = 1:strlength(des_string)
		x_des(i,1) = str2double(des_string(i));
	end
end

