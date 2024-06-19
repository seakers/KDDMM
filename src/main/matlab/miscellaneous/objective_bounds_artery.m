%% Get bounds of objectives for the runs of a particular case for the artery problem
clear
close all
clc

%% Define case parameters
prob_truss = false; % if true -> truss problem, if false -> artery problem
constant_radii = true; % if true -> members have constant radii, if false -> members have discrete radii choices
model = "Truss"; % "Fibre" - fibre model, "Truss" - truss model, "Beam" - beam model

num_runs = 30; % number of runs

% partcoll_bools = [int_pen, AOS, bias_init, ACH];
% nodalprop_bools = [int_pen, AOS, bias_init, ACH];
% orient_bools = [int_pen, AOS, bias_init, ACH];
% inters_bools = [int_pen, AOS, bias_init, ACH];

% Case: Epsilon MOEA
partcoll_bools = [false, false, false, false];
nodalprop_bools = [false, false, false, false];
orient_bools = [false, false, false, false];
inters_bools = [false, false, false, false];

%% Parameter definitions (printable designs)
E = 1.8162e6; % Young's Modulus for polymeric material (example: 1.8162 MPa for SIL material)
sel = 10e-3; % Unit square side length (NOT individual truss length) (in m)
r = 250e-6; % Radius for cross-sectional area of (assumed circular) truss members (in m)
A = pi*(r^2); % Cross-sectional area of truss member
sidenum = 3;
n_members_total = nchoosek(sidenum^2,2); 
c_target = 0.421;

CA_all = get_CA_all(sidenum);
NC = generateNC(sel, sidenum);

n_members_repeated = 2*nchoosek(sidenum,2);
n_variables = n_members_total - n_members_repeated;

%% Read designs for each run and compute first objective and artery problem constraints in the deviation objective
obj1_bounds_allruns = zeros(num_runs,2); % bounds for each run -> [max, min]
obj2_c22c11_bounds_allruns = zeros(num_runs,2);
obj2_c12c11_bounds_allruns = zeros(num_runs,2);
obj2_c21c11_bounds_allruns = zeros(num_runs,2);
obj2_c61_bounds_allruns = zeros(num_runs,2);
obj2_c62_bounds_allruns = zeros(num_runs,2);
obj2_c16_bounds_allruns = zeros(num_runs,2);
obj2_c26_bounds_allruns = zeros(num_runs,2);
obj2_c66c11_bounds_allruns = zeros(num_runs,2);

for i = 1:num_runs
	fieldname = strcat('run',num2str(i));
	[data_array_run, design_array_run] = read_csv_data(prob_truss, model, constant_radii, partcoll_bools, nodalprop_bools, orient_bools, inters_bools, n_members_total, i-1);
	n_des = size(design_array_run,1);
	obj1_vals_run = zeros(n_des, 1);
	obj2_c22c11_vals_run = zeros(n_des, 1);
	obj2_c12c11_vals_run = zeros(n_des, 1);
	obj2_c21c11_vals_run = zeros(n_des, 1);
	obj2_c61_vals_run = zeros(n_des, 1);
	obj2_c62_vals_run = zeros(n_des, 1);
	obj2_c16_vals_run = zeros(n_des, 1);
	obj2_c26_vals_run = zeros(n_des, 1);
	obj2_c66c11_vals_run = zeros(n_des, 1);
    des_count = 1;
	
	for j = 1:n_des
		x_des_current = get_binary_array_from_bitstring(design_array_run{j});
		CA_des = CA_all(x_des_current == 1,:);
		rvar_des = r.*ones(1,size(CA_des,1));
        switch model
            case "Fibre"
                if truss_problem
                    [C11,C22,volfrac_des] = fiberStiffnessModel_rVar_V3(sel,rvar_des,E,CA_des,nucFac,sidenum);
                    C_des = zeros(6);
                    C_des(1,1) = C11;
                    C_des(2,2) = C22;
                else
                    disp("Fiber stiffness model not suitable for artery problem")
                    exit
                end
            case "Truss"
                [C_des, volfrac_des] = trussMetaCalc_NxN_1UC_rVar_AVar(sidenum,sel,rvar_des,E,CA_des);
            case "Beam"
                C_des = Beam_2D_NxN_PBC(sel,sidenum,r,E,CA_des);
                volfrac_des = calcVF_NxN_feasOnly(CA_des,r,sel,sidenum);
        end
        if (any(C_des(1:2,1:2) < 1,'all') || any(C_des > E,'all') || any(isnan(C_des),'all'))
            continue
        end
		obj1_vals_run(des_count,1) = C_des(1,1)/volfrac_des;
		obj2_c22c11_vals_run(des_count,1) = abs((C_des(2,2)/C_des(1,1)) - c_target);
		obj2_c12c11_vals_run(des_count,1) = abs((C_des(1,2)/C_des(1,1)) - 0.0745);
		obj2_c21c11_vals_run(des_count,1) = abs((C_des(2,1)/C_des(1,1)) - 0.0745);
		obj2_c61_vals_run(des_count,1) = abs(C_des(3,1));
		obj2_c62_vals_run(des_count,1) = abs(C_des(3,2));
		obj2_c16_vals_run(des_count,1) = abs(C_des(1,3));
		obj2_c26_vals_run(des_count,1) = abs(C_des(2,3));
		obj2_c66c11_vals_run(des_count,1) = abs((C_des(3,3)/C_des(1,1)) - 5.038);
        
        des_count = des_count + 1;
	end
	
	obj1_bounds_allruns(i,:) = [max(obj1_vals_run(1:des_count-1)), min(obj1_vals_run(1:des_count-1))];
	obj2_c12c11_bounds_allruns(i,:) = [max(obj2_c12c11_vals_run(1:des_count-1)), min(obj2_c12c11_vals_run(1:des_count-1))];
	obj2_c21c11_bounds_allruns(i,:) = [max(obj2_c21c11_vals_run(1:des_count-1)), min(obj2_c21c11_vals_run(1:des_count-1))];
	obj2_c22c11_bounds_allruns(i,:) = [max(obj2_c22c11_vals_run(1:des_count-1)), min(obj2_c22c11_vals_run(1:des_count-1))];
	obj2_c61_bounds_allruns(i,:) = [max(obj2_c61_vals_run(1:des_count-1)), min(obj2_c61_vals_run(1:des_count-1))];
	obj2_c62_bounds_allruns(i,:) = [max(obj2_c62_vals_run(1:des_count-1)), min(obj2_c62_vals_run(1:des_count-1))];
	obj2_c16_bounds_allruns(i,:) = [max(obj2_c16_vals_run(1:des_count-1)), min(obj2_c16_vals_run(1:des_count-1))];
	obj2_c26_bounds_allruns(i,:) = [max(obj2_c26_vals_run(1:des_count-1)), min(obj2_c26_vals_run(1:des_count-1))];
	obj2_c66c11_bounds_allruns(i,:) = [max(obj2_c66c11_vals_run(1:des_count-1)), min(obj2_c66c11_vals_run(1:des_count-1))];
end

obj1_bounds = [max(obj1_bounds_allruns(1,:)), min(obj1_bounds_allruns(2,:))];
obj2_c22c11_bounds = [max(obj2_c22c11_bounds_allruns(1,:)), min(obj2_c22c11_bounds_allruns(2,:))];
obj2_c12c11_bounds = [max(obj2_c12c11_bounds_allruns(1,:)), min(obj2_c12c11_bounds_allruns(2,:))];
obj2_c21c11_bounds = [max(obj2_c21c11_bounds_allruns(1,:)), min(obj2_c21c11_bounds_allruns(2,:))];
obj2_c61_bounds = [max(obj2_c61_bounds_allruns(1,:)), min(obj2_c61_bounds_allruns(2,:))];
obj2_c62_bounds = [max(obj2_c62_bounds_allruns(1,:)), min(obj2_c62_bounds_allruns(2,:))];
obj2_c16_bounds = [max(obj2_c16_bounds_allruns(1,:)), min(obj2_c16_bounds_allruns(2,:))];
obj2_c26_bounds = [max(obj2_c26_bounds_allruns(1,:)), min(obj2_c26_bounds_allruns(2,:))];
obj2_c66c11_bounds = [max(obj2_c66c11_bounds_allruns(1,:)), min(obj2_c66c11_bounds_allruns(2,:))];

%% Functions

function [data_array, design_array] = read_csv_data(problem_truss, choice_of_model, constrad_read, partcoll_bools, nodalprop_bools, orient_bools, inters_bools, n_total_members, run_num)
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
    
    switch choice_of_model
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