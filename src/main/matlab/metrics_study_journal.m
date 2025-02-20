%% Metrics study (truss and artery problems) - Truss and Beam Models
clear 
close all
clc

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
biasFactor = 1;
collapsibilityBiasFac = 0.5;
choice_of_model = "Truss"; % "Fibre" -> fibre model, "Truss" -> truss model, "Beam" -> beam model

%minimum_support = 0.005;

%% Cases to consider for GA data
constrad_read = true;
% Case - Epsilon MOEA
truss_problem = true; % true -> truss problem, false -> artery problem
case1_partcoll_bools = [false, false, false, false];
case1_nodalprop_bools = [false, false, false, false];
case1_orient_bools = [false, false, false, false];
case1_inters_bools = [false, false, false, false];

%% Generate random architectures, objectives, constraints and heuristics
n_des = 300; % number of architectures to generate
n_runs = 10; % number of runs to generate "n_des" architectures

n_total_members = nchoosek(sidenum^2,2);
n_repeated_members = 2*nchoosek(sidenum,2);
n_variables = n_total_members - n_repeated_members;

add_ga_data = false;

bool_des_all = zeros(n_total_members,n_des,n_runs);
bool_des_map = containers.Map;

%orient_avg_all = zeros(n_des,n_runs);
orient_avg_norm_all = struct;

min_dist_pen_pf_all = struct;
min_dist_true_pf_all = struct;

obj1_all = struct;
obj2_all = struct;
obj1_true_all = struct;
obj2_true_all = struct;
obj1_pen_all = struct;
obj2_pen_all = struct;

feas_all = struct;
conn_all = struct;
if truss_problem
    stiff_rat_all = struct;
end

coll_all = struct;
nod_all = struct;
%orient_all = zeros(n_des,n_runs);
orient_norm_all = struct;
inters_all = struct;

support_full_feas_runs = zeros(n_runs,1);
support_full_conn_runs = zeros(n_runs,1);
if truss_problem
    support_full_stiffrat_runs = zeros(n_runs,1);
end
support_full_coll_runs = zeros(n_runs,1);
support_full_nod_runs = zeros(n_runs,1);
support_full_orient_runs = zeros(n_runs,1);
support_full_inters_runs = zeros(n_runs,1);

n_pen_pf_runs = zeros(n_runs,1);
%n_true_pf_runs = zeros(n_runs,1);

support_pen_pf_runs = zeros(n_runs,1);
%support_true_pf_runs = zeros(n_runs,1);

if add_ga_data
    pop_size = 100;

    %[f_pen_combined_case1, f_true_combined_case1, constr_combined_case1, heur_combined_case1, ~] = obtain_combined_data_case(truss_problem, choice_of_model, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_inters_bools, constrad_read, sidenum, pop_size, n_runs);
    
    [f_pen_nonans_allcases, f_true_nonans_allcases, constr_nonans_allcases, heur_nonans_allcases, ~] = obtain_combined_data_allruns(truss_problem, choice_of_model, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_inters_bools, constrad_read, sidenum, pop_size, n_runs);
    
    %f_pen_combined_allcases = f_pen_combined_case1;
    %f_true_combined_allcases = f_true_combined_case1;
    %constr_combined_allcases = constr_combined_case1;
    %heur_combined_allcases = heur_combined_case1;
end

for i = 1:n_runs
    des_count = 1;
    obj1_run = zeros(n_des,1);
    obj2_run = zeros(n_des,1);
    obj1_pen_run = zeros(n_des,1);
    obj2_pen_run = zeros(n_des,1);
    feas_run = zeros(n_des,1);
    if truss_problem
        stiff_rat_run = zeros(n_des,1);
    end
    conn_run = zeros(n_des,1);
    coll_run = zeros(n_des,1);
    nod_run = zeros(n_des,1);
    orient_norm_run = zeros(n_des,1);
    inters_run = zeros(n_des,1);
    bool_des_run = zeros(n_total_members,n_des);
    while des_count<=n_des
        bool_des = randi([0,1],n_variables,1);
        complete_bool_des = get_complete_boolean_array(bool_des, CA_all, sidenum);
        CA_des = CA_all(complete_bool_des~=0,:);
        rvar_des = r.*ones(1,size(CA_des,1));
        switch choice_of_model
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
        valid_design = true;
        if (any(isnan(C_des),'all'))
            valid_design = false;
        end
        if (any(C_des >= E,'all') || any(C_des <= 1,'all'))
            valid_design = false;
        end    
        if ~valid_design
            continue
        else
            if (map_contains_design(bool_des_map,bool_des))
                continue
            else
                bool_des_run(:,des_count) = complete_bool_des;
                bool_des_map(strcat('run',num2str(i),'des',num2str(des_count))) = bool_des_run;
                
                % Constraints
                feas_run(des_count) = feasibility_checker_nonbinary_V5(NC,CA_des,sel,sidenum);
                conn_run(des_count) = connectivityConstraint_PBC_2D(sidenum,NC,CA_des,sel,biasFactor);
                %inters_run(des_count) = feas_run(des_count);
                
                % True Objectives
                if truss_problem
                    obj1_run(des_count) = C_des(2,2);
                    stiff_rat_run(des_count) = abs((C_des(2,2)/C_des(1,1)) - c_ratio);
                    obj2_run(des_count) = volfrac_des;
                else
                    obj1_run(des_count) = C_des(1,1)/volfrac_des;
                    obj2_run(des_count) = (abs((C_des(2,2)/C_des(1,1)) - 0.421) + abs((C_des(2,1)/C_des(1,1)) - 0.0745) + abs((C_des(1,2)/C_des(1,1)) - 0.0745) + abs(C_des(3,1)/E) + abs(C_des(1,3)/E) + abs(C_des(3,2)/E) + abs(C_des(2,3)/E) + abs((C_des(3,3)/C_des(1,1)) - 5.038))/8;
                end
                
                % Penalized Objectives
                obj1_pen_run(des_count) = -obj1_run(des_count)./E;
                obj2_pen_run(des_count) = obj2_run(des_count);
                    
                % Heuristics
                coll_run(des_count) = partCollapseHeuristic_2D(sidenum,CA_des,NC,sel,collapsibilityBiasFac);
                nod_run(des_count) = connectivityHeuristic_2D(sidenum,NC,CA_des,sel,biasFactor);
                %[orient_des, ~] = orientationHeuristic(NC,CA_des,c_ratio);
                [orient_norm_des, ~] = orientationHeuristicNorm(NC,CA_des,sel,c_ratio);
                %orient_run(des_count) = orient_des;
                orient_norm_run(des_count) = orient_norm_des;
                %orient_avg_run(des_count) = avg_angle_des;
                %orient_avg_norm_run(des_count) = avg_angle_norm_des;
                inters_run(des_count) = intersectHeuristic(NC,CA_des);
                
                des_count = des_count + 1;
                disp(strcat(num2str(des_count)," designs generated for run ",num2str(i)))
            end
        end
    end
    disp(strcat("Random designs for run ",num2str(i)," generated"))
    
    if add_ga_data
        current_field = strcat('run',num2str(i));
        f_pen_nonans_currentcase = f_pen_nonans_allcases.(current_field);
        f_true_nonans_currentcase = f_true_nonans_allcases.(current_field);
        constr_nonans_currentcase = constr_nonans_allcases.(current_field);
        heur_nonans_currentcase = heur_nonans_allcases.(current_field);
        
        % Penalized objectives
        obj1_pen_total = cat(1,obj1_pen_run,f_pen_nonans_currentcase(:,1));
        obj2_pen_total = cat(1,obj2_pen_run,f_pen_nonans_currentcase(:,2));
        
        % True objectives
        obj1_total = cat(1,obj1_run,f_true_nonans_currentcase(:,1));
        obj2_total = cat(1,obj2_run,f_true_nonans_currentcase(:,2));
        
        % Constraints
        feas_total = cat(1,feas_run,constr_nonans_currentcase(:,1));
        conn_total = cat(1,conn_run,constr_nonans_currentcase(:,2));
        if truss_problem
            stiffrat_total = cat(1,stiff_rat_run,constr_nonans_currentcase(:,3));
        end
        
        % Heuristics
        coll_total = cat(1,coll_run,heur_nonans_currentcase(:,1));
        nod_total = cat(1,nod_run,heur_nonans_currentcase(:,2));
        orient_total = cat(1,orient_norm_run,heur_nonans_currentcase(:,3));
        %inters_total = feas_total;
        inters_total = cat(1,inters_run,heur_nonans_currentcase(:,4));
    
    else
        % Penalized objectives
        obj1_pen_total = obj1_pen_run;
        obj2_pen_total = obj2_pen_run;
        
        % True objectives
        obj1_total = obj1_run;
        obj2_total = obj2_run;
        
        % Constraints
        feas_total = feas_run;
        conn_total = conn_run;
        if truss_problem
            stiffrat_total = stiff_rat_run;
        end
        
        % Heuristics
        coll_total = coll_run;
        nod_total = nod_run;
        orient_total = orient_norm_run;
        %inters_total = feas_run;
        inters_total = inters_run;
    end
    
    current_field = strcat('trial',num2str(i));
    
    % Constraints penalty values
    feas_all.(current_field) = feas_total;
    conn_all.(current_field) = conn_total;
    
    obj1_true_all.(current_field) = obj1_total;
    obj2_true_all.(current_field) = obj2_total;
    
    pen_objs_pareto = compute_pareto_front(obj1_pen_total,obj2_pen_total);
    true_objs_pareto = compute_pareto_front(-obj1_total,obj2_total);
    true_objs_pareto_correct = [-true_objs_pareto(:,1),true_objs_pareto(:,2)];
    
    %%%% Normalizing objectives and pfs wrt max and min from objectives
    obj1_max = max(obj1_total);
    obj1_min = min(obj1_total);
    obj1_pen_max = max(obj1_pen_total);
    obj1_pen_min = min(obj1_pen_total);
    
    obj2_max = max(obj2_total);
    obj2_min = min(obj2_total);
    obj2_pen_max = max(obj2_pen_total);
    obj2_pen_min = min(obj2_pen_total);
    
    obj1_pen_norm_total = (obj1_pen_total - obj1_pen_min)/(obj1_pen_max - obj1_pen_min);
    obj2_pen_norm_total = (obj2_pen_total - obj2_pen_min)/(obj2_pen_max - obj2_pen_min);
    
    obj1_norm_total = (obj1_total - obj1_min)/(obj1_max - obj1_min);
    obj2_norm_total = (obj2_total - obj2_min)/(obj2_max - obj2_min);
    
    obj1_pen_all.(current_field) = obj1_pen_norm_total;
    obj2_pen_all.(current_field) = obj2_pen_norm_total;
    
    obj1_all.(current_field) = obj1_norm_total;
    obj2_all.(current_field) = obj2_norm_total;
    
    pen_objs_norm_pareto = [(pen_objs_pareto(:,1) - obj1_pen_min)/(obj1_pen_max - obj1_pen_min), (pen_objs_pareto(:,2) - obj2_pen_min)/(obj2_pen_max - obj2_pen_min)]; 
    true_objs_norm_pareto = [(true_objs_pareto_correct(:,1) - obj1_min)/(obj1_max - obj1_min), (true_objs_pareto_correct(:,2) - obj2_min)/(obj2_max - obj2_min)]; 
    
    min_dist_pen_pf_total = zeros(size(obj1_pen_total,1),1);
    min_dist_true_pf_total = zeros(size(obj1_pen_total,1),1);
    for k = 1:size(obj1_pen_total,1)
        min_dist_pen_pf_total(k,1) = compute_min_pf_dist([obj1_pen_norm_total(k,1),obj2_pen_norm_total(k,1)],pen_objs_norm_pareto);
        min_dist_true_pf_total(k,1) = compute_min_pf_dist([obj1_norm_total(k,1),obj2_norm_total(k,1)],true_objs_norm_pareto);
    end
    
    bool_des_all(:,:,i) = bool_des_run;
    
    % Min. dist. to true/penallized Pareto Fronts
    min_dist_pen_pf_all.(current_field) = min_dist_pen_pf_total;
    min_dist_true_pf_all.(current_field) = min_dist_true_pf_total;
    
    % Heuristics
    coll_all.(current_field) = coll_total;
    nod_all.(current_field) = nod_total;
    %orient_all(:,i) = orient_run;
    orient_norm_all.(current_field) = orient_total;
    inters_all.(current_field) = inters_total;
    %orient_avg_all(:,i) = orient_avg_run;
    %orient_avg_norm_all(:,i) = orient_avg_norm_run;
    
    %%% Support calculation
    % Constraints
    support_full_feas_runs(i,1) = length(feas_total(feas_total==1))/size(obj1_total,1);
    support_full_conn_runs(i,1) = length(conn_total(conn_total==1))/size(obj1_total,1);
    if truss_problem
        support_full_stiffrat_runs(i,1) = length(stiffrat_total(stiffrat_total==0))/size(obj1_total,1);
    end
    
    % Heuristics
    support_full_coll_runs(i,1) = length(coll_total(coll_total==1))/size(obj1_total,1);
    support_full_nod_runs(i,1) = length(nod_total(nod_total==1))/size(obj1_total,1);
    support_full_orient_runs(i,1) = length(orient_total(orient_total==1))/size(obj1_total,1);
    support_full_inters_runs(i,1) = length(inters_total(inters_total==1))/size(obj1_total,1);
    
    n_pen_pf_runs(i,1) = size(pen_objs_pareto,1);
    
    support_pen_pf_runs(i,1) = size(pen_objs_pareto,1)/size(obj1_total,1);
end

if truss_problem
    support_tablestats = [mean(support_pen_pf_runs), std(support_pen_pf_runs);
        mean(support_full_feas_runs), std(support_full_feas_runs);
        mean(support_full_conn_runs), std(support_full_conn_runs);
        mean(support_full_stiffrat_runs), std(support_full_stiffrat_runs);
        mean(support_full_coll_runs), std(support_full_coll_runs);
        mean(support_full_nod_runs), std(support_full_nod_runs);
        mean(support_full_orient_runs), std(support_full_orient_runs);
        mean(support_full_inters_runs), std(support_full_inters_runs)];
else
    support_tablestats = [mean(support_pen_pf_runs), std(support_pen_pf_runs);
        mean(support_full_feas_runs), std(support_full_feas_runs);
        mean(support_full_conn_runs), std(support_full_conn_runs);
        mean(support_full_coll_runs), std(support_full_coll_runs);
        mean(support_full_nod_runs), std(support_full_nod_runs);
        mean(support_full_orient_runs), std(support_full_orient_runs);
        mean(support_full_inters_runs), std(support_full_inters_runs)];
end

%% Heuristic violation plots

n_des_total = 0;
for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    obj1_total = obj1_all.(current_field);
    n_des_total = n_des_total + size(obj1_total,1);
end

% Objectives
obj1_array = zeros(n_des_total,1);
obj2_array = zeros(n_des_total,1);

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

% Min. dist. to true/penalized Pareto Fronts
min_dist_pen_pf_array = zeros(n_des_total,1);
min_dist_true_pf_array = zeros(n_des_total,1);

index = 1;
for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    % Objectives
    obj1_total = obj1_all.(current_field);
    obj2_total = obj2_all.(current_field);
    
    obj1_true_total = obj1_true_all.(current_field);
    obj2_true_total = obj2_true_all.(current_field);
    
    % Constraints
    feas_total = feas_all.(current_field);
    conn_total = conn_all.(current_field);
    if truss_problem
        stiffrat_total = stiff_rat_all.(current_field);
    end
    
    % Heuristics
    coll_total = coll_all.(current_field);
    nod_total = nod_all.(current_field);
    orient_total = orient_norm_all.(current_field);
    inters_total = inters_all.(current_field);
    
    % Min. dist. to true/penalized Pareto Fronts
    min_dist_pen_pf_total = min_dist_pen_pf_all.(current_field);
    min_dist_true_pf_total = min_dist_true_pf_all.(current_field);
    
    n_des_run = size(obj1_total,1);
    
    % Populate Objectives
    obj1_array(index:index+n_des_run-1,1) = obj1_total;
    obj2_array(index:index+n_des_run-1,1) = obj2_total;
    
    obj1_true_array(index:index+n_des_run-1,1) = obj1_true_total;
    obj2_true_array(index:index+n_des_run-1,1) = obj2_true_total;
    
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
    
    % Populate Min. dist. to true/penalized Pareto Fronts
    min_dist_pen_pf_array(index:index+n_des_run-1,1) = min_dist_pen_pf_total;
    min_dist_true_pf_array(index:index+n_des_run-1,1) = min_dist_true_pf_total;

    index = index + n_des_run;
end

%%% Determine threshold values for association rule mining
% Objectives
obj1_thresh_val = prctile(obj1_array,75);
obj2_thresh_val = prctile(obj2_array,25);

% Constraints
feas_thresh_val = prctile(feas_array,75);
conn_thresh_val = prctile(conn_array,75);
if truss_problem
    stiffrat_thresh_val = prctile(stiffrat_array,25);
end

% Heuristics
coll_thresh_val = prctile(coll_array,75);
nod_thresh_val = prctile(nod_array,75);
orient_thresh_val = prctile(orient_array,75);
inters_thresh_val = prctile(inters_array,75);

% Min. dist. to true/penalized Pareto Fronts
mindist_pfpen_thresh_val = prctile(min_dist_pen_pf_array,60);
mindist_pftrue_thresh_val = prctile(min_dist_true_pf_array,70);

%%% Computing penalized objectives
% Consttaints
if truss_problem
    stiffrat_constr_array = abs(stiffrat_array)/10;
end
feas_pen_array = log10(abs(feas_array))/16;
conn_pen_array = log10(abs(conn_array))/16;

% Penalized objectives
if truss_problem
    obj1_pen_array = -obj1_true_array./E;  
    obj2_pen_array = obj2_true_array;  
else
    obj1_pen_array = -obj1_true_array./E;  
    obj2_pen_array = obj2_true_array;  
end

if truss_problem
    score_tabelstats = [mean(feas_array), std(feas_array);
        mean(conn_array), std(conn_array);
        mean(stiffrat_array), std(stiffrat_array);
        mean(coll_array), std(coll_array);
        mean(nod_array), std(nod_array);
        mean(orient_array), std(orient_array);
        mean(inters_array), std(inters_array)];
else
    score_tabelstats = [mean(feas_array), std(feas_array);
        mean(conn_array), std(conn_array);
        mean(coll_array), std(coll_array);
        mean(nod_array), std(nod_array);
        mean(orient_array), std(orient_array);
        mean(inters_array), std(inters_array)];
end

figure 
scatter(obj1_pen_array,obj2_pen_array,[],1 - coll_array,'filled')
if truss_problem
    xlabel('Penalized $C_{22}$','Interpreter','Latex','FontSize',16)
    ylabel('Penalized $v_f$','Interpreter','Latex','FontSize',16)
else 
    xlabel('Penalized $\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
    ylabel('Penalized deviation','Interpreter','Latex','FontSize',16)
end
colorbar
title('Partial Collapsibility Violation','FontSize',16)

figure 
scatter(obj1_pen_array,obj2_pen_array,[],1 - nod_array,'filled')
if truss_problem
    xlabel('Penalized $C_{22}$','Interpreter','Latex','FontSize',16)
    ylabel('Penalized $v_f$','Interpreter','Latex','FontSize',16)
else
    xlabel('Penalized $\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
    ylabel('Penalized deviation','Interpreter','Latex','FontSize',16)
end
colorbar;
title('Nodal Properties Violation','FontSize',16)

figure 
scatter(obj1_pen_array,obj2_pen_array,[],1 - orient_array,'filled')
if truss_problem
    xlabel('Penalized $C_{22}$','Interpreter','Latex','FontSize',16)
    ylabel('Penalized $v_f$','Interpreter','Latex','FontSize',16)
else
    xlabel('Penalized $\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
    ylabel('Penalized deviation','Interpreter','Latex','FontSize',16)
end
colorbar;
title('Orientation Violation','FontSize',16)

figure 
scatter(obj1_pen_array,obj2_pen_array,[],1 - inters_array,'filled')
if truss_problem
    xlabel('Penalized $C_{22}$','Interpreter','Latex','FontSize',16)
    ylabel('Penalized $v_f$','Interpreter','Latex','FontSize',16)
else
    xlabel('Penalized $\frac{C_{11}}{v_f}$','Interpreter','Latex','FontSize',16)
    ylabel('Penalized deviation','Interpreter','Latex','FontSize',16)
end
colorbar;
title('Intersection Violation','FontSize',16)

figure 
scatter(obj1_true_array,obj2_true_array,[],1 - coll_array,'filled')
if truss_problem
    xlabel('$C_{22}$','Interpreter','Latex')
    ylabel('$v_f$','Interpreter','Latex')
else
    xlabel('$\frac{C_{11}}{v_f}$','Interpreter','Latex')
    ylabel('deviation','Interpreter','Latex')
end
colorbar
title('Partial Collapsibility Violation')

figure 
scatter(obj1_true_array,obj2_true_array,[],1 - nod_array,'filled')
if truss_problem
    xlabel('$C_{22}$','Interpreter','Latex')
    ylabel('$v_f$','Interpreter','Latex')
else
    xlabel('$\frac{C_{11}}{v_f}$','Interpreter','Latex')
    ylabel('deviation','Interpreter','Latex')
end
colorbar;
title('Nodal Properties Violation')

figure 
scatter(obj1_true_array,obj2_true_array,[],1 - orient_array,'filled')
if truss_problem
    xlabel('$C_{22}$','Interpreter','Latex')
    ylabel('$v_f$','Interpreter','Latex')
else
    xlabel('$\frac{C_{11}}{v_f}$','Interpreter','Latex')
    ylabel('deviation','Interpreter','Latex')
end
colorbar;
title('Orientation Violation')

figure 
scatter(obj1_true_array,obj2_true_array,[],1 - inters_array,'filled')
if truss_problem
    xlabel('$C_{22}$','Interpreter','Latex')
    ylabel('$v_f$','Interpreter','Latex')
else
    xlabel('$\frac{C_{11}}{v_f}$','Interpreter','Latex')
    ylabel('deviation','Interpreter','Latex')
end
colorbar;
title('Intersection Violation')

%% Compute correlation coefficients of each heuristic with objectives and constraints
pearson_coll_truepfdist = zeros(n_runs,1);
pearson_coll_penpfdist = zeros(n_runs,1);
pearson_coll_feas = zeros(n_runs,1);
pearson_coll_conn = zeros(n_runs,1);
% pval_coll_truepfdist = zeros(n_runs,1);
% pval_coll_penpfdist = zeros(n_runs,1);
% pval_coll_feas = zeros(n_runs,1);
% pval_coll_conn = zeros(n_runs,1);

pearson_nod_truepfdist = zeros(n_runs,1);
pearson_nod_penpfdist = zeros(n_runs,1);
pearson_nod_feas = zeros(n_runs,1);
pearson_nod_conn = zeros(n_runs,1);
% pval_nod_truepfdist = zeros(n_runs,1);
% pval_nod_penpfdist = zeros(n_runs,1);
% pval_nod_feas = zeros(n_runs,1);
% pval_nod_conn = zeros(n_runs,1);

pearson_orient_truepfdist = zeros(n_runs,1);
pearson_orient_penpfdist = zeros(n_runs,1);
pearson_orient_feas = zeros(n_runs,1);
pearson_orient_conn = zeros(n_runs,1);
% pval_orient_truepfdist = zeros(n_runs,1);
% pval_orient_penpfdist = zeros(n_runs,1);
% pval_orient_feas = zeros(n_runs,1);
% pval_orient_conn = zeros(n_runs,1);

pearson_inters_truepfdist = zeros(n_runs,1);
pearson_inters_penpfdist = zeros(n_runs,1);
pearson_inters_feas = zeros(n_runs,1);
pearson_inters_conn = zeros(n_runs,1);
% pval_inters_truepfdist = zeros(n_runs,1);
% pval_inters_penpfdist = zeros(n_runs,1);
% pval_inters_feas = zeros(n_runs,1);
% pval_inters_conn = zeros(n_runs,1);

spearman_coll_truepfdist = zeros(n_runs,1);
spearman_coll_penpfdist = zeros(n_runs,1);
spearman_coll_feas = zeros(n_runs,1);
spearman_coll_conn = zeros(n_runs,1);
% pval_coll_truepfdist = zeros(n_runs,1);
% pval_coll_penpfdist = zeros(n_runs,1);
% pval_coll_feas = zeros(n_runs,1);
% pval_coll_conn = zeros(n_runs,1);

spearman_nod_truepfdist = zeros(n_runs,1);
spearman_nod_penpfdist = zeros(n_runs,1);
spearman_nod_feas = zeros(n_runs,1);
spearman_nod_conn = zeros(n_runs,1);
% pval_nod_truepfdist = zeros(n_runs,1);
% pval_nod_penpfdist = zeros(n_runs,1);
% pval_nod_feas = zeros(n_runs,1);
% pval_nod_conn = zeros(n_runs,1);

spearman_orient_truepfdist = zeros(n_runs,1);
spearman_orient_penpfdist = zeros(n_runs,1);
spearman_orient_feas = zeros(n_runs,1);
spearman_orient_conn = zeros(n_runs,1);
% pval_orient_truepfdist = zeros(n_runs,1);
% pval_orient_penpfdist = zeros(n_runs,1);
% pval_orient_feas = zeros(n_runs,1);
% pval_orient_conn = zeros(n_runs,1);

spearman_inters_truepfdist = zeros(n_runs,1);
spearman_inters_penpfdist = zeros(n_runs,1);
spearman_inters_feas = zeros(n_runs,1);
spearman_inters_conn = zeros(n_runs,1);
% pval_inters_truepfdist = zeros(n_runs,1);
% pval_inters_penpfdist = zeros(n_runs,1);
% pval_inters_feas = zeros(n_runs,1);
% pval_inters_conn = zeros(n_runs,1);

if truss_problem
    pearson_coll_stiffrat = zeros(n_runs,1);
    pearson_nod_stiffrat = zeros(n_runs,1);
    pearson_orient_stiffrat = zeros(n_runs,1);
    pearson_inters_stiffrat = zeros(n_runs,1);
    spearman_coll_stiffrat = zeros(n_runs,1);
    spearman_nod_stiffrat = zeros(n_runs,1);
    spearman_orient_stiffrat = zeros(n_runs,1);
    spearman_inters_stiffrat = zeros(n_runs,1);
    % pval_coll_stiffrat = zeros(n_runs,1);
    % pval_nod_stiffrat = zeros(n_runs,1);
    % pval_orient_stiffrat = zeros(n_runs,1);
    % pval_inters_stiffrat = zeros(n_runs,1);
    % pval_coll_stiffrat = zeros(n_runs,1);
    % pval_nod_stiffrat = zeros(n_runs,1);
    % pval_orient_stiffrat = zeros(n_runs,1);
    % pval_inters_stiffrat = zeros(n_runs,1);
end
    
for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    % Min. dist. to true/penalized Pareto Fronts
    min_dist_true_pf_total = min_dist_true_pf_all.(current_field);
    min_dist_pen_pf_total = min_dist_pen_pf_all.(current_field);
    
    % Heuristics
    coll_total = coll_all.(current_field);
    nod_total = nod_all.(current_field);
    orient_total = orient_norm_all.(current_field);
    inters_total = inters_all.(current_field);
    
    % Constraints
    feas_total = feas_all.(current_field);
    conn_total = conn_all.(current_field);
    
    % Computing Pearson's Coefficients
    [pearson_coll_truepfdist(i),~] = corr(coll_total,min_dist_true_pf_total,'Type','Pearson','Rows','complete');
    [pearson_coll_penpfdist(i),~] = corr(coll_total,min_dist_pen_pf_total,'Type','Pearson','Rows','complete');
    [pearson_coll_feas(i),~] = corr(coll_total,feas_total,'Type','Pearson','Rows','complete');
    [pearson_coll_conn(i),~] = corr(coll_total,conn_total,'Type','Pearson','Rows','complete');
    
    [pearson_nod_truepfdist(i),~] = corr(nod_total,min_dist_true_pf_total,'Type','Pearson','Rows','complete');
    [pearson_nod_penpfdist(i),~] = corr(nod_total,min_dist_pen_pf_total,'Type','Pearson','Rows','complete');
    [pearson_nod_feas(i),~] = corr(nod_total,feas_total,'Type','Pearson','Rows','complete');    
    [pearson_nod_conn(i),~] = corr(nod_total,conn_total,'Type','Pearson','Rows','complete');
    
    [pearson_orient_truepfdist(i),~] = corr(orient_total,min_dist_true_pf_total,'Type','Pearson','Rows','complete');
    [pearson_orient_penpfdist(i),~] = corr(orient_total,min_dist_pen_pf_total,'Type','Pearson','Rows','complete');
    [pearson_orient_feas(i),~] = corr(orient_total,feas_total,'Type','Pearson','Rows','complete');   
    [pearson_orient_conn(i),~] = corr(orient_total,conn_total,'Type','Pearson','Rows','complete');
    
    [pearson_inters_truepfdist(i),~] = corr(inters_total,min_dist_true_pf_total,'Type','Pearson','Rows','complete');
    [pearson_inters_penpfdist(i),~] = corr(inters_total,min_dist_pen_pf_total,'Type','Pearson','Rows','complete');
    [pearson_inters_feas(i),~] = corr(inters_total,feas_total,'Type','Pearson','Rows','complete');    
    [pearson_inters_conn(i),~] = corr(inters_total,conn_total,'Type','Pearson','Rows','complete');
    
    % Computing Spearman's Coefficient
    [spearman_coll_truepfdist(i),~] = corr(coll_total,min_dist_true_pf_total,'Type','Spearman','Rows','complete');
    [spearman_coll_penpfdist(i),~] = corr(coll_total,min_dist_pen_pf_total,'Type','Spearman','Rows','complete');
    [spearman_coll_feas(i),~] = corr(coll_total,feas_total,'Type','Spearman','Rows','complete');    
    [spearman_coll_conn(i),~] = corr(coll_total,conn_total,'Type','Spearman','Rows','complete');
    
    [spearman_nod_truepfdist(i),~] = corr(nod_total,min_dist_true_pf_total,'Type','Spearman','Rows','complete');
    [spearman_nod_penpfdist(i),~] = corr(nod_total,min_dist_pen_pf_total,'Type','Spearman','Rows','complete');
    [spearman_nod_feas(i),~] = corr(nod_total,feas_total,'Type','Spearman','Rows','complete');    
    [spearman_nod_conn(i),~] = corr(nod_total,conn_total,'Type','Spearman','Rows','complete');
    
    [spearman_orient_truepfdist(i),~] = corr(orient_total,min_dist_true_pf_total,'Type','Spearman','Rows','complete');
    [spearman_orient_penpfdist(i),~] = corr(orient_total,min_dist_pen_pf_total,'Type','Spearman','Rows','complete');
    [spearman_orient_feas(i),~] = corr(orient_total,feas_total,'Type','Spearman','Rows','complete');    
    [spearman_orient_conn(i),~] = corr(orient_total,conn_total,'Type','Spearman','Rows','complete');
    
    [spearman_inters_truepfdist(i),~] = corr(inters_total,min_dist_true_pf_total,'Type','Spearman','Rows','complete');
    [spearman_inters_penpfdist(i),~] = corr(inters_total,min_dist_pen_pf_total,'Type','Spearman','Rows','complete');
    [spearman_inters_feas(i),~] = corr(inters_total,feas_total,'Type','Spearman','Rows','complete');
    [spearman_inters_conn(i),~] = corr(inters_total,conn_total,'Type','Spearman','Rows','complete');
    
    if truss_problem
        stiffrat_total = stiff_rat_all.(current_field);
        [pearson_coll_stiffrat(i),~] = corr((coll_total),stiffrat_total,'Type','Pearson','Rows','complete');
        [pearson_nod_stiffrat(i),~] = corr(nod_total,stiffrat_total,'Type','Pearson','Rows','complete');
        [pearson_orient_stiffrat(i),~] = corr(orient_total,stiffrat_total,'Type','Pearson','Rows','complete');
        [pearson_inters_stiffrat(i),~] = corr(inters_total,stiffrat_total,'Type','Pearson','Rows','complete');
        [spearman_coll_stiffrat(i),~] = corr((coll_total),stiffrat_total,'Type','Spearman','Rows','complete');
        [spearman_nod_stiffrat(i),~] = corr(nod_total,stiffrat_total,'Type','Spearman','Rows','complete');
        [spearman_orient_stiffrat(i),~] = corr(orient_total,stiffrat_total,'Type','Spearman','Rows','complete');
        [spearman_inters_stiffrat(i),~] = corr(inters_total,stiffrat_total,'Type','Spearman','Rows','complete');
    end
        
end

if truss_problem
    correlation_tablestats = [mean(pearson_coll_truepfdist), std(pearson_coll_truepfdist), mean(spearman_coll_truepfdist), std(spearman_coll_truepfdist);
        mean(pearson_coll_feas), std(pearson_coll_feas), mean(spearman_coll_feas), std(spearman_coll_feas);
        mean(pearson_coll_conn), std(pearson_coll_conn), mean(spearman_coll_conn), std(spearman_coll_conn);
        mean(pearson_coll_stiffrat), std(pearson_coll_stiffrat), mean(spearman_coll_stiffrat), std(spearman_coll_stiffrat);
        mean(pearson_nod_truepfdist), std(pearson_nod_truepfdist), mean(spearman_nod_truepfdist), std(spearman_nod_truepfdist);
        mean(pearson_nod_feas), std(pearson_nod_feas), mean(spearman_nod_feas), std(spearman_nod_feas);
        mean(pearson_nod_conn), std(pearson_nod_conn), mean(spearman_nod_conn), std(spearman_nod_conn);
        mean(pearson_nod_stiffrat), std(pearson_nod_stiffrat), mean(spearman_nod_stiffrat), std(spearman_nod_stiffrat);
        mean(pearson_orient_truepfdist), std(pearson_orient_truepfdist), mean(spearman_orient_truepfdist), std(spearman_orient_truepfdist);
        mean(pearson_orient_feas), std(pearson_orient_feas), mean(spearman_orient_feas), std(spearman_orient_feas);
        mean(pearson_orient_conn), std(pearson_orient_conn), mean(spearman_orient_conn), std(spearman_orient_conn);
        mean(pearson_orient_stiffrat), std(pearson_orient_stiffrat), mean(spearman_orient_stiffrat), std(spearman_orient_stiffrat);
        mean(pearson_inters_truepfdist), std(pearson_inters_truepfdist), mean(spearman_inters_truepfdist), std(spearman_inters_truepfdist);
        mean(pearson_inters_feas), std(pearson_inters_feas), mean(spearman_inters_feas), std(spearman_inters_feas);
        mean(pearson_inters_conn), std(pearson_inters_conn), mean(spearman_inters_conn), std(spearman_inters_conn);
        mean(pearson_inters_stiffrat), std(pearson_inters_stiffrat), mean(spearman_inters_stiffrat), std(spearman_inters_stiffrat)];
else
    correlation_tablestats = [mean(pearson_coll_truepfdist), std(pearson_coll_truepfdist), mean(spearman_coll_truepfdist), std(spearman_coll_truepfdist);
        mean(pearson_coll_feas), std(pearson_coll_feas), mean(spearman_coll_feas), std(spearman_coll_feas);
        mean(pearson_coll_conn), std(pearson_coll_conn), mean(spearman_coll_conn), std(spearman_coll_conn);
        mean(pearson_nod_truepfdist), std(pearson_nod_truepfdist), mean(spearman_nod_truepfdist), std(spearman_nod_truepfdist);
        mean(pearson_nod_feas), std(pearson_nod_feas), mean(spearman_nod_feas), std(spearman_nod_feas);
        mean(pearson_nod_conn), std(pearson_nod_conn), mean(spearman_nod_conn), std(spearman_nod_conn);
        mean(pearson_orient_truepfdist), std(pearson_orient_truepfdist), mean(spearman_orient_truepfdist), std(spearman_orient_truepfdist);
        mean(pearson_orient_feas), std(pearson_orient_feas), mean(spearman_orient_feas), std(spearman_orient_feas);
        mean(pearson_orient_conn), std(pearson_orient_conn), mean(spearman_orient_conn), std(spearman_orient_conn);
        mean(pearson_inters_truepfdist), std(pearson_inters_truepfdist), mean(spearman_inters_truepfdist), std(spearman_inters_truepfdist);
        mean(pearson_inters_feas), std(pearson_inters_feas), mean(spearman_inters_feas), std(spearman_inters_feas);
        mean(pearson_inters_conn), std(pearson_inters_conn), mean(spearman_inters_conn), std(spearman_inters_conn)];
end
   
correlation_penpf_tablestats = [mean(pearson_coll_penpfdist), std(pearson_coll_penpfdist), mean(spearman_coll_penpfdist), std(spearman_coll_penpfdist);
    mean(pearson_nod_penpfdist), std(pearson_nod_penpfdist), mean(spearman_nod_penpfdist), std(spearman_nod_penpfdist);
    mean(pearson_orient_penpfdist), std(pearson_orient_penpfdist), mean(spearman_orient_penpfdist), std(spearman_orient_penpfdist);
    mean(pearson_inters_penpfdist), std(pearson_inters_penpfdist), mean(spearman_inters_penpfdist), std(spearman_inters_penpfdist)];

%% Thresholding heuristics, objectives and constraints into high and low
% Heuristics
coll_all_thresh = struct;
nod_all_thresh = struct;
orient_all_thresh = struct;
inters_all_thresh = struct;

% Objectives
obj1_all_thresh = struct;
obj2_all_thresh = struct;

% Constraints
feas_all_thresh = struct;
conn_all_thresh = struct;
if truss_problem
    stiffrat_all_thresh = struct;
end

% Min. dist. to true/penalized Pareto Fronts
min_dist_pen_pf_all_thresh = struct;
min_dist_true_pf_all_thresh = struct;

for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    % Heuristics
    coll_total = coll_all.(current_field);
    nod_total = nod_all.(current_field);
    orient_total = orient_norm_all.(current_field);
    inters_total = inters_all.(current_field);
    
    coll_run_thresh = double(coll_total >= coll_thresh_val);
    nod_run_thresh = double(nod_total >= nod_thresh_val);
    orient_run_thresh = double(orient_total >= orient_thresh_val);
    inters_run_thresh = double(inters_total >= inters_thresh_val);
    
    coll_all_thresh.(current_field) = coll_run_thresh;
    nod_all_thresh.(current_field) = nod_run_thresh;
    orient_all_thresh.(current_field) = orient_run_thresh;
    inters_all_thresh.(current_field) = inters_run_thresh;
    
    % Objectives
    obj1_total = obj1_all.(current_field);
    obj2_total = obj2_all.(current_field);
    
    obj1_run_thresh = double(obj1_total >= obj1_thresh_val);
    obj2_run_thresh = double(obj2_total >= obj2_thresh_val);
    
    obj1_all_thresh.(current_field) = obj1_run_thresh;
    obj2_all_thresh.(current_field) = obj2_run_thresh;
    
    % Constraints
    feas_total = feas_all.(current_field);
    conn_total = conn_all.(current_field);
    
    feas_run_thresh = double(feas_total >= feas_thresh_val);
    conn_run_thresh = double(conn_total >= conn_thresh_val);

    feas_all_thresh.(current_field) = feas_run_thresh;
    conn_all_thresh.(current_field) = conn_run_thresh;
    
    if truss_problem
        stiffrat_total = stiff_rat_all.(current_field);
        stiffrat_run_thresh = double(stiffrat_total >= stiffrat_thresh_val);
        stiffrat_all_thresh.(current_field) = stiffrat_run_thresh;
    end
        
    % Min. dist. to true/penalized Pareto Fronts
    min_dist_pen_pf_total = min_dist_pen_pf_all.(current_field);
    min_dist_true_pf_total = min_dist_true_pf_all.(current_field);

    min_dist_pen_pf_thresh = double(min_dist_pen_pf_total >= mindist_pfpen_thresh_val);
    min_dist_true_pf_thresh = double(min_dist_true_pf_total >= mindist_pftrue_thresh_val);
    
    min_dist_pen_pf_all_thresh.(current_field) = min_dist_pen_pf_thresh;
    min_dist_true_pf_all_thresh.(current_field) = min_dist_true_pf_thresh;
     
end
    
%% Computing interestingness measures for different association rules
% intrness = [support_a, support_b, support_ab, confidence_a2b,
% confidence_b2a, lift]
% heur_intrness_run = [low_heur_intrness; high_heur_intrness]{size = 10 x 6 x 2} one 3D array for each run  

coll_intrness_truepf_constr = struct;
nod_intrness_truepf_constr = struct;
orient_intrness_truepf_constr = struct;
inters_intrness_truepf_constr = struct;

coll_intrness_truepf = struct;
nod_intrness_truepf = struct;
orient_intrness_truepf = struct;
inters_intrness_truepf = struct;

coll_intrness_penpf = struct;
nod_intrness_penpf = struct;
orient_intrness_penpf = struct;
inters_intrness_penpf = struct;

for i = 1:n_runs
    
    current_field = strcat('trial',num2str(i));
    % Objectives
    obj1_thresh = obj1_all_thresh.(current_field);
    obj2_thresh = obj2_all_thresh.(current_field);
    
    % Heuristics
    coll_thresh = coll_all_thresh.(current_field);
    nod_thresh = nod_all_thresh.(current_field);
    orient_thresh = orient_all_thresh.(current_field);
    inters_thresh = inters_all_thresh.(current_field);
    
    % Constraints
    feas_thresh = feas_all_thresh.(current_field);
    conn_thresh = conn_all_thresh.(current_field);
    if truss_problem
        stiffrat_thresh = stiffrat_all_thresh.(current_field);
    end
    
    % Min. dist. to true/penalized Pareto Fronts
    min_dist_pen_pf_thresh = min_dist_pen_pf_all_thresh.(current_field);
    min_dist_true_pf_thresh = min_dist_true_pf_all_thresh.(current_field);
    
    % Compute interestingness measures wrt true Pareto Front distance
    coll_intrness_truepf_run = compute_heur_intrness_penpf(coll_thresh, min_dist_true_pf_thresh);
    coll_intrness_truepf.(current_field) = coll_intrness_truepf_run;
    nod_intrness_truepf_run = compute_heur_intrness_penpf(nod_thresh, min_dist_true_pf_thresh);
    nod_intrness_truepf.(current_field) = nod_intrness_truepf_run;
    orient_intrness_truepf_run = compute_heur_intrness_penpf(orient_thresh, min_dist_true_pf_thresh);
    orient_intrness_truepf.(current_field) = orient_intrness_truepf_run;
    inters_intrness_truepf_run = compute_heur_intrness_penpf(inters_thresh, min_dist_true_pf_thresh);
    inters_intrness_truepf.(current_field) = inters_intrness_truepf_run;
    
    % Compute interestingness measures wrt to constraints and parts of true
    % Pareto Front
    if truss_problem
        coll_intrness_run = compute_heur_intrness_array_truss(coll_thresh, min_dist_true_pf_thresh, obj1_thresh, obj2_thresh, feas_thresh, conn_thresh, stiffrat_thresh);
        nod_intrness_run = compute_heur_intrness_array_truss(nod_thresh, min_dist_true_pf_thresh, obj1_thresh, obj2_thresh, feas_thresh, conn_thresh, stiffrat_thresh);
        orient_intrness_run = compute_heur_intrness_array_truss(orient_thresh, min_dist_true_pf_thresh, obj1_thresh, obj2_thresh, feas_thresh, conn_thresh, stiffrat_thresh);
        inters_intrness_run = compute_heur_intrness_array_truss(inters_thresh, min_dist_true_pf_thresh, obj1_thresh, obj2_thresh, feas_thresh, conn_thresh, stiffrat_thresh);
    else
        coll_intrness_run = compute_heur_intrness_array_artery(coll_thresh, min_dist_true_pf_thresh, obj1_thresh, obj2_thresh, feas_thresh, conn_thresh);
        nod_intrness_run = compute_heur_intrness_array_artery(nod_thresh, min_dist_true_pf_thresh, obj1_thresh, obj2_thresh, feas_thresh, conn_thresh);
        orient_intrness_run = compute_heur_intrness_array_artery(orient_thresh, min_dist_true_pf_thresh, obj1_thresh, obj2_thresh, feas_thresh, conn_thresh);
        inters_intrness_run = compute_heur_intrness_array_artery(inters_thresh, min_dist_true_pf_thresh, obj1_thresh, obj2_thresh, feas_thresh, conn_thresh);
    end
    
    coll_intrness_truepf_constr.(current_field) = coll_intrness_run;
    nod_intrness_truepf_constr.(current_field) = nod_intrness_run;
    orient_intrness_truepf_constr.(current_field) = orient_intrness_run;
    inters_intrness_truepf_constr.(current_field) = inters_intrness_run;
    
    % Compute interestingnes measures wrt penalized Pareto Front distance
    coll_intrness_penpf_run = compute_heur_intrness_penpf(coll_thresh, min_dist_pen_pf_thresh);
    coll_intrness_penpf.(current_field) = coll_intrness_penpf_run;
    nod_intrness_penpf_run = compute_heur_intrness_penpf(nod_thresh, min_dist_pen_pf_thresh);
    nod_intrness_penpf.(current_field) = nod_intrness_penpf_run;
    orient_intrness_penpf_run = compute_heur_intrness_penpf(orient_thresh, min_dist_pen_pf_thresh);
    orient_intrness_penpf.(current_field) = orient_intrness_penpf_run;
    inters_intrness_penpf_run = compute_heur_intrness_penpf(inters_thresh, min_dist_pen_pf_thresh);
    inters_intrness_penpf.(current_field) = inters_intrness_penpf_run;
    
end
    
%% Obtain interestingness measure arrays for each heuristic in turns
low_heur_close2truepf_lowobj1_runs = zeros(n_runs,6);
low_heur_close2truepf_highobj1_runs = zeros(n_runs,6);
low_heur_close2truepf_lowobj2_runs = zeros(n_runs,6);
low_heur_close2truepf_highobj2_runs = zeros(n_runs,6);
low_heur_low_feas_runs = zeros(n_runs,6);
low_heur_high_feas_runs = zeros(n_runs,6);
low_heur_low_conn_runs = zeros(n_runs,6);
low_heur_high_conn_runs = zeros(n_runs,6);

high_heur_close2truepf_lowobj1_runs = zeros(n_runs,6);
high_heur_close2truepf_highobj1_runs = zeros(n_runs,6);
high_heur_close2truepf_lowobj2_runs = zeros(n_runs,6);
high_heur_close2truepf_highobj2_runs = zeros(n_runs,6);
high_heur_low_feas_runs = zeros(n_runs,6);
high_heur_high_feas_runs = zeros(n_runs,6);
high_heur_low_conn_runs = zeros(n_runs,6);
high_heur_high_conn_runs = zeros(n_runs,6);

if truss_problem
    low_heur_low_stiffrat_runs = zeros(n_runs,6);
    low_heur_high_stiffrat_runs = zeros(n_runs,6);
    high_heur_low_stiffrat_runs = zeros(n_runs,6);
    high_heur_high_stiffrat_runs = zeros(n_runs,6);
end

low_heur_close2penpf = zeros(n_runs,6);
high_heur_close2penpf = zeros(n_runs,6);

low_heur_close2truepf = zeros(n_runs,6);
high_heur_close2truepf = zeros(n_runs,6);

for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    %heur_intrness_run = coll_intrness_truepf_constr.(current_field); % size: 10 x 6 x 2
    %heur_intrness_run = nod_intrness_truepf_constr.(current_field);
    %heur_intrness_run = orient_intrness_truepf_constr.(current_field);
    heur_intrness_run = inters_intrness_truepf_constr.(current_field); 
    
    low_heur_close2truepf_lowobj1_runs(i,:) = heur_intrness_run(1,:,1);
    low_heur_close2truepf_highobj1_runs(i,:) = heur_intrness_run(2,:,1);
    low_heur_close2truepf_lowobj2_runs(i,:) = heur_intrness_run(3,:,1);
    low_heur_close2truepf_highobj2_runs(i,:) = heur_intrness_run(4,:,1);
    low_heur_low_feas_runs(i,:) = heur_intrness_run(5,:,1);
    low_heur_high_feas_runs(i,:) = heur_intrness_run(6,:,1);
    low_heur_low_conn_runs(i,:) = heur_intrness_run(7,:,1);
    low_heur_high_conn_runs(i,:) = heur_intrness_run(8,:,1);
    
    high_heur_close2truepf_lowobj1_runs(i,:) = heur_intrness_run(1,:,2);
    high_heur_close2truepf_highobj1_runs(i,:) = heur_intrness_run(2,:,2);
    high_heur_close2truepf_lowobj2_runs(i,:) = heur_intrness_run(3,:,2);
    high_heur_close2truepf_highobj2_runs(i,:) = heur_intrness_run(4,:,2);
    high_heur_low_feas_runs(i,:) = heur_intrness_run(5,:,2);
    high_heur_high_feas_runs(i,:) = heur_intrness_run(6,:,2);
    high_heur_low_conn_runs(i,:) = heur_intrness_run(7,:,2);
    high_heur_high_conn_runs(i,:) = heur_intrness_run(8,:,2);
    
    if truss_problem
        low_heur_low_stiffrat_runs(i,:) = heur_intrness_run(9,:,1);
        low_heur_high_stiffrat_runs(i,:) = heur_intrness_run(10,:,1);
        high_heur_low_stiffrat_runs(i,:) = heur_intrness_run(9,:,2);
        high_heur_high_stiffrat_runs(i,:) = heur_intrness_run(10,:,2);
    end
    
    %heur_intrness_penpf_run = coll_intrness_penpf.(current_field); % size: 2 x 6
    %heur_intrness_penpf_run = nod_intrness_penpf.(current_field);
    %heur_intrness_penpf_run = orient_intrness_penpf.(current_field);
    heur_intrness_penpf_run = inters_intrness_penpf.(current_field);
    
    %heur_intrness_truepf_run = coll_intrness_truepf.(current_field); % size: 2 x 6
    %heur_intrness_truepf_run = nod_intrness_truepf.(current_field);
    %heur_intrness_truepf_run = orient_intrness_truepf.(current_field);
    heur_intrness_truepf_run = inters_intrness_truepf.(current_field);
    
    low_heur_close2penpf(i,:) = heur_intrness_penpf_run(1,:);
    high_heur_close2penpf(i,:) = heur_intrness_penpf_run(2,:);
    
    low_heur_close2truepf(i,:) = heur_intrness_truepf_run(1,:);
    high_heur_close2truepf(i,:) = heur_intrness_truepf_run(2,:);
    
end

if truss_problem
    low_heur_tablestats = [mean(low_heur_close2truepf_lowobj1_runs(:,1)), std(low_heur_close2truepf_lowobj1_runs(:,1)), mean(low_heur_close2truepf_lowobj1_runs(:,2)), std(low_heur_close2truepf_lowobj1_runs(:,2)), mean(low_heur_close2truepf_lowobj1_runs(:,3)), std(low_heur_close2truepf_lowobj1_runs(:,3)), mean(low_heur_close2truepf_lowobj1_runs(:,4)), std(low_heur_close2truepf_lowobj1_runs(:,4)), mean(low_heur_close2truepf_lowobj1_runs(:,5)), std(low_heur_close2truepf_lowobj1_runs(:,5)), mean(low_heur_close2truepf_lowobj1_runs(:,6)), std(low_heur_close2truepf_lowobj1_runs(:,6));
        mean(low_heur_close2truepf_highobj1_runs(:,1)), std(low_heur_close2truepf_highobj1_runs(:,1)), mean(low_heur_close2truepf_highobj1_runs(:,2)), std(low_heur_close2truepf_highobj1_runs(:,2)), mean(low_heur_close2truepf_highobj1_runs(:,3)), std(low_heur_close2truepf_highobj1_runs(:,3)), mean(low_heur_close2truepf_highobj1_runs(:,4)), std(low_heur_close2truepf_highobj1_runs(:,4)), mean(low_heur_close2truepf_highobj1_runs(:,5)), std(low_heur_close2truepf_highobj1_runs(:,5)), mean(low_heur_close2truepf_highobj1_runs(:,6)), std(low_heur_close2truepf_highobj1_runs(:,6));
        mean(low_heur_close2truepf_lowobj2_runs(:,1)), std(low_heur_close2truepf_lowobj2_runs(:,1)), mean(low_heur_close2truepf_lowobj2_runs(:,2)), std(low_heur_close2truepf_lowobj2_runs(:,2)), mean(low_heur_close2truepf_lowobj2_runs(:,3)), std(low_heur_close2truepf_lowobj2_runs(:,3)), mean(low_heur_close2truepf_lowobj2_runs(:,4)), std(low_heur_close2truepf_lowobj2_runs(:,4)), mean(low_heur_close2truepf_lowobj2_runs(:,5)), std(low_heur_close2truepf_lowobj2_runs(:,5)), mean(low_heur_close2truepf_lowobj2_runs(:,6)), std(low_heur_close2truepf_lowobj2_runs(:,6));
        mean(low_heur_close2truepf_highobj2_runs(:,1)), std(low_heur_close2truepf_highobj2_runs(:,1)), mean(low_heur_close2truepf_highobj2_runs(:,2)), std(low_heur_close2truepf_highobj2_runs(:,2)), mean(low_heur_close2truepf_highobj2_runs(:,3)), std(low_heur_close2truepf_highobj2_runs(:,3)), mean(low_heur_close2truepf_highobj2_runs(:,4)), std(low_heur_close2truepf_highobj2_runs(:,4)), mean(low_heur_close2truepf_highobj2_runs(:,5)), std(low_heur_close2truepf_highobj2_runs(:,5)), mean(low_heur_close2truepf_highobj2_runs(:,6)), std(low_heur_close2truepf_highobj2_runs(:,6));
        mean(low_heur_low_feas_runs(:,1)), std(low_heur_low_feas_runs(:,1)), mean(low_heur_low_feas_runs(:,2)), std(low_heur_low_feas_runs(:,2)), mean(low_heur_low_feas_runs(:,3)), std(low_heur_low_feas_runs(:,3)), mean(low_heur_low_feas_runs(:,4)), std(low_heur_low_feas_runs(:,4)), mean(low_heur_low_feas_runs(:,5)), std(low_heur_low_feas_runs(:,5)), mean(low_heur_low_feas_runs(:,6)), std(low_heur_low_feas_runs(:,6));
        mean(low_heur_high_feas_runs(:,1)), std(low_heur_high_feas_runs(:,1)), mean(low_heur_high_feas_runs(:,2)), std(low_heur_high_feas_runs(:,2)), mean(low_heur_high_feas_runs(:,3)), std(low_heur_high_feas_runs(:,3)), mean(low_heur_high_feas_runs(:,4)), std(low_heur_high_feas_runs(:,4)), mean(low_heur_high_feas_runs(:,5)), std(low_heur_high_feas_runs(:,5)), mean(low_heur_high_feas_runs(:,6)), std(low_heur_high_feas_runs(:,6));
        mean(low_heur_low_conn_runs(:,1)), std(low_heur_low_conn_runs(:,1)), mean(low_heur_low_conn_runs(:,2)), std(low_heur_low_conn_runs(:,2)), mean(low_heur_low_conn_runs(:,3)), std(low_heur_low_conn_runs(:,3)), mean(low_heur_low_conn_runs(:,4)), std(low_heur_low_conn_runs(:,4)), mean(low_heur_low_conn_runs(:,5)), std(low_heur_low_conn_runs(:,5)), mean(low_heur_low_conn_runs(:,6)), std(low_heur_low_conn_runs(:,6));
        mean(low_heur_high_conn_runs(:,1)), std(low_heur_high_conn_runs(:,1)), mean(low_heur_high_conn_runs(:,2)), std(low_heur_high_conn_runs(:,2)), mean(low_heur_high_conn_runs(:,3)), std(low_heur_high_conn_runs(:,3)), mean(low_heur_high_conn_runs(:,4)), std(low_heur_high_conn_runs(:,4)), mean(low_heur_high_conn_runs(:,5)), std(low_heur_high_conn_runs(:,5)), mean(low_heur_high_conn_runs(:,6)), std(low_heur_high_conn_runs(:,6));
        mean(low_heur_low_stiffrat_runs(:,1)), std(low_heur_low_stiffrat_runs(:,1)), mean(low_heur_low_stiffrat_runs(:,2)), std(low_heur_low_stiffrat_runs(:,2)), mean(low_heur_low_stiffrat_runs(:,3)), std(low_heur_low_stiffrat_runs(:,3)), mean(low_heur_low_stiffrat_runs(:,4)), std(low_heur_low_stiffrat_runs(:,4)), mean(low_heur_low_stiffrat_runs(:,5)), std(low_heur_low_stiffrat_runs(:,5)), mean(low_heur_low_stiffrat_runs(:,6)), std(low_heur_low_stiffrat_runs(:,6));
        mean(low_heur_high_stiffrat_runs(:,1)), std(low_heur_high_stiffrat_runs(:,1)), mean(low_heur_high_stiffrat_runs(:,2)), std(low_heur_high_stiffrat_runs(:,2)), mean(low_heur_high_stiffrat_runs(:,3)), std(low_heur_high_stiffrat_runs(:,3)), mean(low_heur_high_stiffrat_runs(:,4)), std(low_heur_high_stiffrat_runs(:,4)), mean(low_heur_high_stiffrat_runs(:,5)), std(low_heur_high_stiffrat_runs(:,5)), mean(low_heur_high_stiffrat_runs(:,6)), std(low_heur_high_stiffrat_runs(:,6))];

    high_heur_tablestats = [mean(high_heur_close2truepf_lowobj1_runs(:,1)), std(high_heur_close2truepf_lowobj1_runs(:,1)), mean(high_heur_close2truepf_lowobj1_runs(:,2)), std(high_heur_close2truepf_lowobj1_runs(:,2)), mean(high_heur_close2truepf_lowobj1_runs(:,3)), std(high_heur_close2truepf_lowobj1_runs(:,3)), mean(high_heur_close2truepf_lowobj1_runs(:,4)), std(high_heur_close2truepf_lowobj1_runs(:,4)), mean(high_heur_close2truepf_lowobj1_runs(:,5)), std(high_heur_close2truepf_lowobj1_runs(:,5)), mean(high_heur_close2truepf_lowobj1_runs(:,6)), std(high_heur_close2truepf_lowobj1_runs(:,6));
        mean(high_heur_close2truepf_highobj1_runs(:,1)), std(high_heur_close2truepf_highobj1_runs(:,1)), mean(high_heur_close2truepf_highobj1_runs(:,2)), std(high_heur_close2truepf_highobj1_runs(:,2)), mean(high_heur_close2truepf_highobj1_runs(:,3)), std(high_heur_close2truepf_highobj1_runs(:,3)), mean(high_heur_close2truepf_highobj1_runs(:,4)), std(high_heur_close2truepf_highobj1_runs(:,4)), mean(high_heur_close2truepf_highobj1_runs(:,5)), std(high_heur_close2truepf_highobj1_runs(:,5)), mean(high_heur_close2truepf_highobj1_runs(:,6)), std(high_heur_close2truepf_highobj1_runs(:,6));
        mean(high_heur_close2truepf_lowobj2_runs(:,1)), std(high_heur_close2truepf_lowobj2_runs(:,1)), mean(high_heur_close2truepf_lowobj2_runs(:,2)), std(high_heur_close2truepf_lowobj2_runs(:,2)), mean(high_heur_close2truepf_lowobj2_runs(:,3)), std(high_heur_close2truepf_lowobj2_runs(:,3)), mean(high_heur_close2truepf_lowobj2_runs(:,4)), std(high_heur_close2truepf_lowobj2_runs(:,4)), mean(high_heur_close2truepf_lowobj2_runs(:,5)), std(high_heur_close2truepf_lowobj2_runs(:,5)), mean(high_heur_close2truepf_lowobj2_runs(:,6)), std(high_heur_close2truepf_lowobj2_runs(:,6));
        mean(high_heur_close2truepf_highobj2_runs(:,1)), std(high_heur_close2truepf_highobj2_runs(:,1)), mean(high_heur_close2truepf_highobj2_runs(:,2)), std(high_heur_close2truepf_highobj2_runs(:,2)), mean(high_heur_close2truepf_highobj2_runs(:,3)), std(high_heur_close2truepf_highobj2_runs(:,3)), mean(high_heur_close2truepf_highobj2_runs(:,4)), std(high_heur_close2truepf_highobj2_runs(:,4)), mean(high_heur_close2truepf_highobj2_runs(:,5)), std(high_heur_close2truepf_highobj2_runs(:,5)), mean(high_heur_close2truepf_highobj2_runs(:,6)), std(high_heur_close2truepf_highobj2_runs(:,6));
        mean(high_heur_low_feas_runs(:,1)), std(high_heur_low_feas_runs(:,1)), mean(high_heur_low_feas_runs(:,2)), std(high_heur_low_feas_runs(:,2)), mean(high_heur_low_feas_runs(:,3)), std(high_heur_low_feas_runs(:,3)), mean(high_heur_low_feas_runs(:,4)), std(high_heur_low_feas_runs(:,4)), mean(high_heur_low_feas_runs(:,5)), std(high_heur_low_feas_runs(:,5)), mean(high_heur_low_feas_runs(:,6)), std(high_heur_low_feas_runs(:,6));
        mean(high_heur_high_feas_runs(:,1)), std(high_heur_high_feas_runs(:,1)), mean(high_heur_high_feas_runs(:,2)), std(high_heur_high_feas_runs(:,2)), mean(high_heur_high_feas_runs(:,3)), std(high_heur_high_feas_runs(:,3)), mean(high_heur_high_feas_runs(:,4)), std(high_heur_high_feas_runs(:,4)), mean(high_heur_high_feas_runs(:,5)), std(high_heur_high_feas_runs(:,5)), mean(high_heur_high_feas_runs(:,6)), std(high_heur_high_feas_runs(:,6));
        mean(high_heur_low_conn_runs(:,1)), std(high_heur_low_conn_runs(:,1)), mean(high_heur_low_conn_runs(:,2)), std(high_heur_low_conn_runs(:,2)), mean(high_heur_low_conn_runs(:,3)), std(high_heur_low_conn_runs(:,3)), mean(high_heur_low_conn_runs(:,4)), std(high_heur_low_conn_runs(:,4)), mean(high_heur_low_conn_runs(:,5)), std(high_heur_low_conn_runs(:,5)), mean(high_heur_low_conn_runs(:,6)), std(high_heur_low_conn_runs(:,6));
        mean(high_heur_high_conn_runs(:,1)), std(high_heur_high_conn_runs(:,1)), mean(high_heur_high_conn_runs(:,2)), std(high_heur_high_conn_runs(:,2)), mean(high_heur_high_conn_runs(:,3)), std(high_heur_high_conn_runs(:,3)), mean(high_heur_high_conn_runs(:,4)), std(high_heur_high_conn_runs(:,4)), mean(high_heur_high_conn_runs(:,5)), std(high_heur_high_conn_runs(:,5)), mean(high_heur_high_conn_runs(:,6)), std(high_heur_high_conn_runs(:,6));
        mean(high_heur_low_stiffrat_runs(:,1)), std(high_heur_low_stiffrat_runs(:,1)), mean(high_heur_low_stiffrat_runs(:,2)), std(high_heur_low_stiffrat_runs(:,2)), mean(high_heur_low_stiffrat_runs(:,3)), std(high_heur_low_stiffrat_runs(:,3)), mean(high_heur_low_stiffrat_runs(:,4)), std(high_heur_low_stiffrat_runs(:,4)), mean(high_heur_low_stiffrat_runs(:,5)), std(high_heur_low_stiffrat_runs(:,5)), mean(high_heur_low_stiffrat_runs(:,6)), std(high_heur_low_stiffrat_runs(:,6));
        mean(high_heur_high_stiffrat_runs(:,1)), std(high_heur_high_stiffrat_runs(:,1)), mean(high_heur_high_stiffrat_runs(:,2)), std(high_heur_high_stiffrat_runs(:,2)), mean(high_heur_high_stiffrat_runs(:,3)), std(high_heur_high_stiffrat_runs(:,3)), mean(high_heur_high_stiffrat_runs(:,4)), std(high_heur_high_stiffrat_runs(:,4)), mean(high_heur_high_stiffrat_runs(:,5)), std(high_heur_high_stiffrat_runs(:,5)), mean(high_heur_high_stiffrat_runs(:,6)), std(high_heur_high_stiffrat_runs(:,6))];
    
else
    low_heur_tablestats = [mean(low_heur_close2truepf_lowobj1_runs(:,1)), std(low_heur_close2truepf_lowobj1_runs(:,1)), mean(low_heur_close2truepf_lowobj1_runs(:,2)), std(low_heur_close2truepf_lowobj1_runs(:,2)), mean(low_heur_close2truepf_lowobj1_runs(:,3)), std(low_heur_close2truepf_lowobj1_runs(:,3)), mean(low_heur_close2truepf_lowobj1_runs(:,4)), std(low_heur_close2truepf_lowobj1_runs(:,4)), mean(low_heur_close2truepf_lowobj1_runs(:,5)), std(low_heur_close2truepf_lowobj1_runs(:,5)), mean(low_heur_close2truepf_lowobj1_runs(:,6)), std(low_heur_close2truepf_lowobj1_runs(:,6));
        mean(low_heur_close2truepf_highobj1_runs(:,1)), std(low_heur_close2truepf_highobj1_runs(:,1)), mean(low_heur_close2truepf_highobj1_runs(:,2)), std(low_heur_close2truepf_highobj1_runs(:,2)), mean(low_heur_close2truepf_highobj1_runs(:,3)), std(low_heur_close2truepf_highobj1_runs(:,3)), mean(low_heur_close2truepf_highobj1_runs(:,4)), std(low_heur_close2truepf_highobj1_runs(:,4)), mean(low_heur_close2truepf_highobj1_runs(:,5)), std(low_heur_close2truepf_highobj1_runs(:,5)), mean(low_heur_close2truepf_highobj1_runs(:,6)), std(low_heur_close2truepf_highobj1_runs(:,6));
        mean(low_heur_close2truepf_lowobj2_runs(:,1)), std(low_heur_close2truepf_lowobj2_runs(:,1)), mean(low_heur_close2truepf_lowobj2_runs(:,2)), std(low_heur_close2truepf_lowobj2_runs(:,2)), mean(low_heur_close2truepf_lowobj2_runs(:,3)), std(low_heur_close2truepf_lowobj2_runs(:,3)), mean(low_heur_close2truepf_lowobj2_runs(:,4)), std(low_heur_close2truepf_lowobj2_runs(:,4)), mean(low_heur_close2truepf_lowobj2_runs(:,5)), std(low_heur_close2truepf_lowobj2_runs(:,5)), mean(low_heur_close2truepf_lowobj2_runs(:,6)), std(low_heur_close2truepf_lowobj2_runs(:,6));
        mean(low_heur_close2truepf_highobj2_runs(:,1)), std(low_heur_close2truepf_highobj2_runs(:,1)), mean(low_heur_close2truepf_highobj2_runs(:,2)), std(low_heur_close2truepf_highobj2_runs(:,2)), mean(low_heur_close2truepf_highobj2_runs(:,3)), std(low_heur_close2truepf_highobj2_runs(:,3)), mean(low_heur_close2truepf_highobj2_runs(:,4)), std(low_heur_close2truepf_highobj2_runs(:,4)), mean(low_heur_close2truepf_highobj2_runs(:,5)), std(low_heur_close2truepf_highobj2_runs(:,5)), mean(low_heur_close2truepf_highobj2_runs(:,6)), std(low_heur_close2truepf_highobj2_runs(:,6));
        mean(low_heur_low_feas_runs(:,1)), std(low_heur_low_feas_runs(:,1)), mean(low_heur_low_feas_runs(:,2)), std(low_heur_low_feas_runs(:,2)), mean(low_heur_low_feas_runs(:,3)), std(low_heur_low_feas_runs(:,3)), mean(low_heur_low_feas_runs(:,4)), std(low_heur_low_feas_runs(:,4)), mean(low_heur_low_feas_runs(:,5)), std(low_heur_low_feas_runs(:,5)), mean(low_heur_low_feas_runs(:,6)), std(low_heur_low_feas_runs(:,6));
        mean(low_heur_high_feas_runs(:,1)), std(low_heur_high_feas_runs(:,1)), mean(low_heur_high_feas_runs(:,2)), std(low_heur_high_feas_runs(:,2)), mean(low_heur_high_feas_runs(:,3)), std(low_heur_high_feas_runs(:,3)), mean(low_heur_high_feas_runs(:,4)), std(low_heur_high_feas_runs(:,4)), mean(low_heur_high_feas_runs(:,5)), std(low_heur_high_feas_runs(:,5)), mean(low_heur_high_feas_runs(:,6)), std(low_heur_high_feas_runs(:,6));
        mean(low_heur_low_conn_runs(:,1)), std(low_heur_low_conn_runs(:,1)), mean(low_heur_low_conn_runs(:,2)), std(low_heur_low_conn_runs(:,2)), mean(low_heur_low_conn_runs(:,3)), std(low_heur_low_conn_runs(:,3)), mean(low_heur_low_conn_runs(:,4)), std(low_heur_low_conn_runs(:,4)), mean(low_heur_low_conn_runs(:,5)), std(low_heur_low_conn_runs(:,5)), mean(low_heur_low_conn_runs(:,6)), std(low_heur_low_conn_runs(:,6));
        mean(low_heur_high_conn_runs(:,1)), std(low_heur_high_conn_runs(:,1)), mean(low_heur_high_conn_runs(:,2)), std(low_heur_high_conn_runs(:,2)), mean(low_heur_high_conn_runs(:,3)), std(low_heur_high_conn_runs(:,3)), mean(low_heur_high_conn_runs(:,4)), std(low_heur_high_conn_runs(:,4)), mean(low_heur_high_conn_runs(:,5)), std(low_heur_high_conn_runs(:,5)), mean(low_heur_high_conn_runs(:,6)), std(low_heur_high_conn_runs(:,6))];

    high_heur_tablestats = [mean(high_heur_close2truepf_lowobj1_runs(:,1)), std(high_heur_close2truepf_lowobj1_runs(:,1)), mean(high_heur_close2truepf_lowobj1_runs(:,2)), std(high_heur_close2truepf_lowobj1_runs(:,2)), mean(high_heur_close2truepf_lowobj1_runs(:,3)), std(high_heur_close2truepf_lowobj1_runs(:,3)), mean(high_heur_close2truepf_lowobj1_runs(:,4)), std(high_heur_close2truepf_lowobj1_runs(:,4)), mean(high_heur_close2truepf_lowobj1_runs(:,5)), std(high_heur_close2truepf_lowobj1_runs(:,5)), mean(high_heur_close2truepf_lowobj1_runs(:,6)), std(high_heur_close2truepf_lowobj1_runs(:,6));
        mean(high_heur_close2truepf_highobj1_runs(:,1)), std(high_heur_close2truepf_highobj1_runs(:,1)), mean(high_heur_close2truepf_highobj1_runs(:,2)), std(high_heur_close2truepf_highobj1_runs(:,2)), mean(high_heur_close2truepf_highobj1_runs(:,3)), std(high_heur_close2truepf_highobj1_runs(:,3)), mean(high_heur_close2truepf_highobj1_runs(:,4)), std(high_heur_close2truepf_highobj1_runs(:,4)), mean(high_heur_close2truepf_highobj1_runs(:,5)), std(high_heur_close2truepf_highobj1_runs(:,5)), mean(high_heur_close2truepf_highobj1_runs(:,6)), std(high_heur_close2truepf_highobj1_runs(:,6));
        mean(high_heur_close2truepf_lowobj2_runs(:,1)), std(high_heur_close2truepf_lowobj2_runs(:,1)), mean(high_heur_close2truepf_lowobj2_runs(:,2)), std(high_heur_close2truepf_lowobj2_runs(:,2)), mean(high_heur_close2truepf_lowobj2_runs(:,3)), std(high_heur_close2truepf_lowobj2_runs(:,3)), mean(high_heur_close2truepf_lowobj2_runs(:,4)), std(high_heur_close2truepf_lowobj2_runs(:,4)), mean(high_heur_close2truepf_lowobj2_runs(:,5)), std(high_heur_close2truepf_lowobj2_runs(:,5)), mean(high_heur_close2truepf_lowobj2_runs(:,6)), std(high_heur_close2truepf_lowobj2_runs(:,6));
        mean(high_heur_close2truepf_highobj2_runs(:,1)), std(high_heur_close2truepf_highobj2_runs(:,1)), mean(high_heur_close2truepf_highobj2_runs(:,2)), std(high_heur_close2truepf_highobj2_runs(:,2)), mean(high_heur_close2truepf_highobj2_runs(:,3)), std(high_heur_close2truepf_highobj2_runs(:,3)), mean(high_heur_close2truepf_highobj2_runs(:,4)), std(high_heur_close2truepf_highobj2_runs(:,4)), mean(high_heur_close2truepf_highobj2_runs(:,5)), std(high_heur_close2truepf_highobj2_runs(:,5)), mean(high_heur_close2truepf_highobj2_runs(:,6)), std(high_heur_close2truepf_highobj2_runs(:,6));
        mean(high_heur_low_feas_runs(:,1)), std(high_heur_low_feas_runs(:,1)), mean(high_heur_low_feas_runs(:,2)), std(high_heur_low_feas_runs(:,2)), mean(high_heur_low_feas_runs(:,3)), std(high_heur_low_feas_runs(:,3)), mean(high_heur_low_feas_runs(:,4)), std(high_heur_low_feas_runs(:,4)), mean(high_heur_low_feas_runs(:,5)), std(high_heur_low_feas_runs(:,5)), mean(high_heur_low_feas_runs(:,6)), std(high_heur_low_feas_runs(:,6));
        mean(high_heur_high_feas_runs(:,1)), std(high_heur_high_feas_runs(:,1)), mean(high_heur_high_feas_runs(:,2)), std(high_heur_high_feas_runs(:,2)), mean(high_heur_high_feas_runs(:,3)), std(high_heur_high_feas_runs(:,3)), mean(high_heur_high_feas_runs(:,4)), std(high_heur_high_feas_runs(:,4)), mean(high_heur_high_feas_runs(:,5)), std(high_heur_high_feas_runs(:,5)), mean(high_heur_high_feas_runs(:,6)), std(high_heur_high_feas_runs(:,6));
        mean(high_heur_low_conn_runs(:,1)), std(high_heur_low_conn_runs(:,1)), mean(high_heur_low_conn_runs(:,2)), std(high_heur_low_conn_runs(:,2)), mean(high_heur_low_conn_runs(:,3)), std(high_heur_low_conn_runs(:,3)), mean(high_heur_low_conn_runs(:,4)), std(high_heur_low_conn_runs(:,4)), mean(high_heur_low_conn_runs(:,5)), std(high_heur_low_conn_runs(:,5)), mean(high_heur_low_conn_runs(:,6)), std(high_heur_low_conn_runs(:,6));
        mean(high_heur_high_conn_runs(:,1)), std(high_heur_high_conn_runs(:,1)), mean(high_heur_high_conn_runs(:,2)), std(high_heur_high_conn_runs(:,2)), mean(high_heur_high_conn_runs(:,3)), std(high_heur_high_conn_runs(:,3)), mean(high_heur_high_conn_runs(:,4)), std(high_heur_high_conn_runs(:,4)), mean(high_heur_high_conn_runs(:,5)), std(high_heur_high_conn_runs(:,5)), mean(high_heur_high_conn_runs(:,6)), std(high_heur_high_conn_runs(:,6))];
end

heur_penpf_tablestats = [mean(low_heur_close2penpf(:,1)), std(low_heur_close2penpf(:,1)), mean(low_heur_close2penpf(:,2)), std(low_heur_close2penpf(:,2)), mean(low_heur_close2penpf(:,3)), std(low_heur_close2penpf(:,3)), mean(low_heur_close2penpf(:,4)), std(low_heur_close2penpf(:,4)), mean(low_heur_close2penpf(:,5)), std(low_heur_close2penpf(:,5)), mean(low_heur_close2penpf(:,6)), std(low_heur_close2penpf(:,6));
    mean(high_heur_close2penpf(:,1)), std(high_heur_close2penpf(:,1)), mean(high_heur_close2penpf(:,2)), std(high_heur_close2penpf(:,2)), mean(high_heur_close2penpf(:,3)), std(high_heur_close2penpf(:,3)), mean(high_heur_close2penpf(:,4)), std(high_heur_close2penpf(:,4)), mean(high_heur_close2penpf(:,5)), std(high_heur_close2penpf(:,5)), mean(high_heur_close2penpf(:,6)), std(high_heur_close2penpf(:,6))];

heur_truepf_tablestats = [mean(low_heur_close2truepf(:,1)), std(low_heur_close2truepf(:,1)), mean(low_heur_close2truepf(:,2)), std(low_heur_close2truepf(:,2)), mean(low_heur_close2truepf(:,3)), std(low_heur_close2truepf(:,3)), mean(low_heur_close2truepf(:,4)), std(low_heur_close2truepf(:,4)), mean(low_heur_close2truepf(:,5)), std(low_heur_close2truepf(:,5)), mean(low_heur_close2truepf(:,6)), std(low_heur_close2truepf(:,6));
    mean(high_heur_close2truepf(:,1)), std(high_heur_close2truepf(:,1)), mean(high_heur_close2truepf(:,2)), std(high_heur_close2truepf(:,2)), mean(high_heur_close2truepf(:,3)), std(high_heur_close2truepf(:,3)), mean(high_heur_close2truepf(:,4)), std(high_heur_close2truepf(:,4)), mean(high_heur_close2truepf(:,5)), std(high_heur_close2truepf(:,5)), mean(high_heur_close2truepf(:,6)), std(high_heur_close2truepf(:,6))];

%% Functions - CHANGE 
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

function plot_regression(linmodel,X,y,Xname,Yname,run_num)
    X_norm = (X - mean(X))/std(X);
    Y_pred_norm = predict(linmodel,X_norm);
    Y_pred = (Y_pred_norm*std(y)) + mean(y);
    figure
    plot(X, y, '*b')
    hold on
    plot(X, Y_pred, 'k')
    hold off
    xlabel(Xname,'Interpreter','Latex')
    ylabel(Yname,'Interpreter','Latex')
    title(strcat('Linear Regression fit for run ',num2str(run_num)))
end

function contains = map_contains_design(design_map, design) 
    contains = false;
    for j = keys(design_map)
        key = j{1};
        if isequal(design_map(key),design)
            contains = true;
            break
        end
    end
end

function heur_intrness_array = compute_heur_intrness_array_truss(heur_thresh, mindist_truepf_thresh, c22thresh, volfracthresh, feasthresh, connthresh, stiffratthresh)
    heur_intrness_array = zeros(10,6,2);
    
    [sup_lowheur_close2truepf_lowc22, conf_lowheur_close2truepf_lowc22, lift_lowheur_close2truepf_lowc22] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_truepf_thresh & ~c22thresh), sum(~heur_thresh & ~mindist_truepf_thresh & ~c22thresh), size(c22thresh,1));
    [sup_lowheur_close2truepf_highc22, conf_lowheur_close2truepf_highc22, lift_lowheur_close2truepf_highc22] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_truepf_thresh & c22thresh), sum(~heur_thresh & ~mindist_truepf_thresh & c22thresh), size(c22thresh,1));
    [sup_lowheur_close2truepf_lowvf, conf_lowheur_close2truepf_lowvf, lift_lowheur_close2truepf_lowvf] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_truepf_thresh & ~volfracthresh), sum(~heur_thresh & ~mindist_truepf_thresh & ~volfracthresh), size(c22thresh,1));
    [sup_lowheur_close2truepf_highvf, conf_lowcoll_close2truepf_highvf, lift_lowheur_close2truepf_highvf] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_truepf_thresh & volfracthresh), sum(~heur_thresh & ~mindist_truepf_thresh & volfracthresh), size(c22thresh,1));
    [sup_lowheur_lowfeas, conf_lowheur_lowfeas, lift_lowheur_lowfeas] = compute_interestingness_measures(sum(~heur_thresh), sum(~feasthresh), sum(~heur_thresh & ~feasthresh), size(c22thresh,1));
    [sup_lowheur_highfeas, conf_lowheur_highfeas, lift_lowheur_highfeas] = compute_interestingness_measures(sum(~heur_thresh), sum(feasthresh), sum(~heur_thresh & feasthresh), size(c22thresh,1));
    [sup_lowheur_lowconn, conf_lowheur_lowconn, lift_lowheur_lowconn] = compute_interestingness_measures(sum(~heur_thresh), sum(~connthresh), sum(~heur_thresh & ~connthresh), size(c22thresh,1));
    [sup_lowheur_highconn, conf_lowheur_highconn, lift_lowheur_highconn] = compute_interestingness_measures(sum(~heur_thresh), sum(connthresh), sum(~heur_thresh & connthresh), size(c22thresh,1));
    [sup_lowheur_lowstiffrat, conf_lowheur_lowstiffrat, lift_lowheur_lowstiffrat] = compute_interestingness_measures(sum(~heur_thresh), sum(~stiffratthresh), sum(~heur_thresh & ~stiffratthresh), size(c22thresh,1));
    [sup_lowheur_highstiffrat, conf_lowheur_highstiffrat, lift_lowheur_highstiffrat] = compute_interestingness_measures(sum(~heur_thresh), sum(stiffratthresh), sum(~heur_thresh & stiffratthresh), size(c22thresh,1));
    
    [sup_highheur_close2truepf_lowc22, conf_highheur_close2truepf_lowc22, lift_highheur_close2truepf_lowc22] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_truepf_thresh & ~c22thresh), sum(heur_thresh & ~mindist_truepf_thresh & ~c22thresh), size(c22thresh,1));
    [sup_highheur_close2truepf_highc22, conf_highheur_close2truepf_highc22, lift_highheur_close2truepf_highc22] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_truepf_thresh & c22thresh), sum(heur_thresh & ~mindist_truepf_thresh & c22thresh), size(c22thresh,1));
    [sup_highheur_close2truepf_lowvf, conf_highheur_close2truepf_lowvf, lift_highheur_close2truepf_lowvf] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_truepf_thresh & ~volfracthresh), sum(heur_thresh & ~mindist_truepf_thresh & ~volfracthresh), size(c22thresh,1));
    [sup_highheur_close2truepf_highvf, conf_highheur_close2truepf_highvf, lift_highheur_close2truepf_highvf] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_truepf_thresh & volfracthresh), sum(heur_thresh & ~mindist_truepf_thresh & volfracthresh), size(c22thresh,1));
    [sup_highheur_lowfeas, conf_highheur_lowfeas, lift_highheur_lowfeas] = compute_interestingness_measures(sum(heur_thresh), sum(~feasthresh), sum(heur_thresh & ~feasthresh), size(c22thresh,1));
    [sup_highheur_highfeas, conf_highheur_highfeas, lift_highheur_highfeas] = compute_interestingness_measures(sum(heur_thresh), sum(feasthresh), sum(heur_thresh & feasthresh), size(c22thresh,1));
    [sup_highheur_lowconn, conf_highheur_lowconn, lift_highheur_lowconn] = compute_interestingness_measures(sum(heur_thresh), sum(~connthresh), sum(heur_thresh & ~connthresh), size(c22thresh,1));
    [sup_highheur_highconn, conf_highheur_highconn, lift_highheur_highconn] = compute_interestingness_measures(sum(heur_thresh), sum(connthresh), sum(heur_thresh & connthresh), size(c22thresh,1));
    [sup_highheur_lowstiffrat, conf_highheur_lowstiffrat, lift_highheur_lowstiffrat] = compute_interestingness_measures(sum(heur_thresh), sum(~stiffratthresh), sum(heur_thresh & ~stiffratthresh), size(c22thresh,1));
    [sup_highheur_highstiffrat, conf_highheur_highstiffrat, lift_highheur_highstiffrat] = compute_interestingness_measures(sum(heur_thresh), sum(stiffratthresh), sum(heur_thresh & stiffratthresh), size(c22thresh,1));
    
    heur_intrness_array(:,:,1) = [sup_lowheur_close2truepf_lowc22, conf_lowheur_close2truepf_lowc22, lift_lowheur_close2truepf_lowc22;
                            sup_lowheur_close2truepf_highc22, conf_lowheur_close2truepf_highc22, lift_lowheur_close2truepf_highc22;
                            sup_lowheur_close2truepf_lowvf, conf_lowheur_close2truepf_lowvf, lift_lowheur_close2truepf_lowvf;
                            sup_lowheur_close2truepf_highvf, conf_lowcoll_close2truepf_highvf, lift_lowheur_close2truepf_highvf;
                            sup_lowheur_lowfeas, conf_lowheur_lowfeas, lift_lowheur_lowfeas;
                            sup_lowheur_highfeas, conf_lowheur_highfeas, lift_lowheur_highfeas;
                            sup_lowheur_lowconn, conf_lowheur_lowconn, lift_lowheur_lowconn;
                            sup_lowheur_highconn, conf_lowheur_highconn, lift_lowheur_highconn;
                            sup_lowheur_lowstiffrat, conf_lowheur_lowstiffrat, lift_lowheur_lowstiffrat;
                            sup_lowheur_highstiffrat, conf_lowheur_highstiffrat, lift_lowheur_highstiffrat];
                        
     heur_intrness_array(:,:,2) = [sup_highheur_close2truepf_lowc22, conf_highheur_close2truepf_lowc22, lift_highheur_close2truepf_lowc22;
                            sup_highheur_close2truepf_highc22, conf_highheur_close2truepf_highc22, lift_highheur_close2truepf_highc22;
                            sup_highheur_close2truepf_lowvf, conf_highheur_close2truepf_lowvf, lift_highheur_close2truepf_lowvf;
                            sup_highheur_close2truepf_highvf, conf_highheur_close2truepf_highvf, lift_highheur_close2truepf_highvf;
                            sup_highheur_lowfeas, conf_highheur_lowfeas, lift_highheur_lowfeas;
                            sup_highheur_highfeas, conf_highheur_highfeas, lift_highheur_highfeas;
                            sup_highheur_lowconn, conf_highheur_lowconn, lift_highheur_lowconn;
                            sup_highheur_highconn, conf_highheur_highconn, lift_highheur_highconn;
                            sup_highheur_lowstiffrat, conf_highheur_lowstiffrat, lift_highheur_lowstiffrat;
                            sup_highheur_highstiffrat, conf_highheur_highstiffrat, lift_highheur_highstiffrat];
                        
end

function heur_intrness_array = compute_heur_intrness_array_artery(heur_thresh, mindist_truepf_thresh, stiffthresh, devthresh, feasthresh, connthresh)
    heur_intrness_array = zeros(8,6,2);
    
    [sup_lowheur_close2truepf_lowstiff, conf_lowheur_close2truepf_lowstiff, lift_lowheur_close2truepf_lowstiff] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_truepf_thresh & ~stiffthresh), sum(~heur_thresh & ~mindist_truepf_thresh & ~stiffthresh), size(stiffthresh,1));
    [sup_lowheur_close2truepf_highstiff, conf_lowheur_close2truepf_highstiff, lift_lowheur_close2truepf_highstiff] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_truepf_thresh & stiffthresh), sum(~heur_thresh & ~mindist_truepf_thresh & stiffthresh), size(stiffthresh,1));
    [sup_lowheur_close2truepf_lowdev, conf_lowheur_close2truepf_lowdev, lift_lowheur_close2truepf_lowdev] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_truepf_thresh & ~devthresh), sum(~heur_thresh & ~mindist_truepf_thresh & ~devthresh), size(stiffthresh,1));
    [sup_lowheur_close2truepf_highdev, conf_lowcoll_close2truepf_highdev, lift_lowheur_close2truepf_highdev] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_truepf_thresh & devthresh), sum(~heur_thresh & ~mindist_truepf_thresh & devthresh), size(stiffthresh,1));
    [sup_lowheur_lowfeas, conf_lowheur_lowfeas, lift_lowheur_lowfeas] = compute_interestingness_measures(sum(~heur_thresh), sum(~feasthresh), sum(~heur_thresh & ~feasthresh), size(stiffthresh,1));
    [sup_lowheur_highfeas, conf_lowheur_highfeas, lift_lowheur_highfeas] = compute_interestingness_measures(sum(~heur_thresh), sum(feasthresh), sum(~heur_thresh & feasthresh), size(stiffthresh,1));
    [sup_lowheur_lowconn, conf_lowheur_lowconn, lift_lowheur_lowconn] = compute_interestingness_measures(sum(~heur_thresh), sum(~connthresh), sum(~heur_thresh & ~connthresh), size(stiffthresh,1));
    [sup_lowheur_highconn, conf_lowheur_highconn, lift_lowheur_highconn] = compute_interestingness_measures(sum(~heur_thresh), sum(connthresh), sum(~heur_thresh & connthresh), size(stiffthresh,1));
    
    [sup_highheur_close2truepf_lowstiff, conf_highheur_close2truepf_lowstiff, lift_highheur_close2truepf_lowstiff] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_truepf_thresh & ~stiffthresh), sum(heur_thresh & ~mindist_truepf_thresh & ~stiffthresh), size(stiffthresh,1));
    [sup_highheur_close2truepf_highstiff, conf_highheur_close2truepf_highstiff, lift_highheur_close2truepf_highstiff] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_truepf_thresh & stiffthresh), sum(heur_thresh & ~mindist_truepf_thresh & stiffthresh), size(stiffthresh,1));
    [sup_highheur_close2truepf_lowdev, conf_highheur_close2truepf_lowdev, lift_highheur_close2truepf_lowdev] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_truepf_thresh & ~devthresh), sum(heur_thresh & ~mindist_truepf_thresh & ~devthresh), size(stiffthresh,1));
    [sup_highheur_close2truepf_highdev, conf_highheur_close2truepf_highdev, lift_highheur_close2truepf_highdev] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_truepf_thresh & devthresh), sum(heur_thresh & ~mindist_truepf_thresh & devthresh), size(stiffthresh,1));
    [sup_highheur_lowfeas, conf_highheur_lowfeas, lift_highheur_lowfeas] = compute_interestingness_measures(sum(heur_thresh), sum(~feasthresh), sum(heur_thresh & ~feasthresh), size(stiffthresh,1));
    [sup_highheur_highfeas, conf_highheur_highfeas, lift_highheur_highfeas] = compute_interestingness_measures(sum(heur_thresh), sum(feasthresh), sum(heur_thresh & feasthresh), size(stiffthresh,1));
    [sup_highheur_lowconn, conf_highheur_lowconn, lift_highheur_lowconn] = compute_interestingness_measures(sum(heur_thresh), sum(~connthresh), sum(heur_thresh & ~connthresh), size(stiffthresh,1));
    [sup_highheur_highconn, conf_highheur_highconn, lift_highheur_highconn] = compute_interestingness_measures(sum(heur_thresh), sum(connthresh), sum(heur_thresh & connthresh), size(stiffthresh,1));
    
    heur_intrness_array(:,:,1) = [sup_lowheur_close2truepf_lowstiff, conf_lowheur_close2truepf_lowstiff, lift_lowheur_close2truepf_lowstiff;
                            sup_lowheur_close2truepf_highstiff, conf_lowheur_close2truepf_highstiff, lift_lowheur_close2truepf_highstiff;
                            sup_lowheur_close2truepf_lowdev, conf_lowheur_close2truepf_lowdev, lift_lowheur_close2truepf_lowdev;
                            sup_lowheur_close2truepf_highdev, conf_lowcoll_close2truepf_highdev, lift_lowheur_close2truepf_highdev;
                            sup_lowheur_lowfeas, conf_lowheur_lowfeas, lift_lowheur_lowfeas;
                            sup_lowheur_highfeas, conf_lowheur_highfeas, lift_lowheur_highfeas;
                            sup_lowheur_lowconn, conf_lowheur_lowconn, lift_lowheur_lowconn;
                            sup_lowheur_highconn, conf_lowheur_highconn, lift_lowheur_highconn];
                        
     heur_intrness_array(:,:,2) = [sup_highheur_close2truepf_lowstiff, conf_highheur_close2truepf_lowstiff, lift_highheur_close2truepf_lowstiff;
                            sup_highheur_close2truepf_highstiff, conf_highheur_close2truepf_highstiff, lift_highheur_close2truepf_highstiff;
                            sup_highheur_close2truepf_lowdev, conf_highheur_close2truepf_lowdev, lift_highheur_close2truepf_lowdev;
                            sup_highheur_close2truepf_highdev, conf_highheur_close2truepf_highdev, lift_highheur_close2truepf_highdev;
                            sup_highheur_lowfeas, conf_highheur_lowfeas, lift_highheur_lowfeas;
                            sup_highheur_highfeas, conf_highheur_highfeas, lift_highheur_highfeas;
                            sup_highheur_lowconn, conf_highheur_lowconn, lift_highheur_lowconn;
                            sup_highheur_highconn, conf_highheur_highconn, lift_highheur_highconn];
                        
end

function heur_intr_array = compute_heur_intrness_penpf(heur_thresh, mindist_penpf_thresh)
    heur_intr_array = zeros(2,6);
    [sup_lowheur_close2penpf, conf_lowheur_close2penpf, lift_lowheur_close2penpf] = compute_interestingness_measures(sum(~heur_thresh), sum(~mindist_penpf_thresh), sum(~heur_thresh & ~mindist_penpf_thresh), size(mindist_penpf_thresh,1));
    [sup_highheur_close2penpf, conf_highheur_close2penpf, lift_highheur_close2penpf] = compute_interestingness_measures(sum(heur_thresh), sum(~mindist_penpf_thresh), sum(heur_thresh & ~mindist_penpf_thresh), size(mindist_penpf_thresh,1));
    
    heur_intr_array(1,:) = [sup_lowheur_close2penpf, conf_lowheur_close2penpf, lift_lowheur_close2penpf];
    heur_intr_array(2,:) = [sup_highheur_close2penpf, conf_highheur_close2penpf, lift_highheur_close2penpf];
end

function [support_vals, confidence_vals, lift_val] = compute_interestingness_measures(n_a,n_b,n_ab,n_all)
    % compute support, confidence and lift vals for association rule A->B
    sup_a = compute_support_arm(n_a,n_all);
    sup_b = compute_support_arm(n_b,n_all);
    sup_ab = compute_support_arm(n_ab,n_all);
    conf_a2b = compute_confidence_arm(n_ab,n_a);
    conf_b2a = compute_confidence_arm(n_ab,n_b);
    lift_val = compute_lift_arm(n_a,n_b,n_ab,n_all);
    support_vals = [sup_a, sup_b, sup_ab];
    confidence_vals = [conf_a2b, conf_b2a];
    
end

function lift = compute_lift_arm(n_X,n_Y,n_XY,n_total) 
    lift = (n_XY/n_X)/(n_Y/n_total);
end

function support = compute_support_arm(n_X,n_total)
    support = n_X/n_total;
end

function confidence = compute_confidence_arm(n_XY,n_X)
    confidence = n_XY/n_X;
end

function targetAngle = compute_target_orientation(target) 
    if target < 1
        ntarget = 1/target;
    else
        ntarget = target;
    end
    k = 5;
    tA = 90./(1+exp(-k.*(ntarget-1)));
    if target < 1
        targetAngle = 90 - tA;
    else
        targetAngle = tA;
    end
end

function complete_bool_array = get_complete_boolean_array(bool_des,CA_full,sidenum)
    top_nodes = get_top_edge_nodes(sidenum);
    %n_total_members = nchoosek(sidenum^2,2);
    complete_bool_array = [];
    n_count = 0;
    for i=1:(sidenum^2)
        for j=i+1:(sidenum^2)
            right_edge = false;
            top_edge = false;
            if (i > (sidenum^2 - sidenum))
                if (j > (sidenum^2 - sidenum))
                    repeated_member = [i - ((sidenum-1)*(sidenum)), j - ((sidenum-1)*sidenum)];
                    repeat_index = find_member_index(CA_full, repeated_member);
                    complete_bool_array = [complete_bool_array;bool_des(repeat_index)];
                    right_edge = true;
                end
            end
            if ismember(i,top_nodes)
                if ismember(j,top_nodes)
                    repeated_member = [i - (sidenum-1), j - (sidenum-1)];
                    repeat_index = find_member_index(CA_full, repeated_member);
                    complete_bool_array = [complete_bool_array;bool_des(repeat_index)];
                    top_edge = true;
                end
            end
            if (~right_edge && ~top_edge)
                complete_bool_array = [complete_bool_array;bool_des(n_count+1)];
                n_count = n_count + 1;
            end
        end
    end        
end

function [objs_pen_nonans_allcases, objs_true_nonans_allcases, constraints_nonans_allcases, heuristics_nonans_allcases, designs_nonans_allcases] = obtain_combined_data_allruns(truss_prob, model_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, constrad_prob_read, sidenodenum, n_pop, n_runs)
    
    n_members_total = nchoosek(sidenodenum^2,2);    
   
    objs_pen_nonans_allcases = struct;
    objs_true_nonans_allcases = struct;
    constraints_nonans_allcases = struct;
    heuristics_nonans_allcases = struct;
    designs_nonans_allcases = struct;
    
    for i = 1:n_runs
        [data_array, designs_array] = read_csv_data(truss_prob, model_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, constrad_prob_read, n_members_total, i-1);
        
        n_constr = size(data_array,2) - 4 - 4; % number of constraints changes based on problem, so subtract 4 (2 pen. and 2 true objs.) and 4 (heurs) from number of total columns
        data_array_nonans_bool = any(isnan(data_array),2);
        data_array_nonans = data_array(~data_array_nonans_bool,:);
        designs_array_nonans = designs_array(~data_array_nonans_bool,:);
        
        current_field = strcat('run',num2str(i));
        
        objs_pen_nonans_allcases.(current_field) = [data_array_nonans(:,1), data_array_nonans(:,2)];
        objs_true_nonans_allcases.(current_field) = [data_array_nonans(:,3), data_array_nonans(:,4)];
        constraints_nonans_allcases.(current_field) = data_array_nonans(:,5:5+n_constr-1);
        heuristics_nonans_allcases.(current_field) = data_array_nonans(:,5+n_constr:end);
        designs_nonans_allcases.(current_field) = designs_array_nonans;
        
    end
    
end

function [objs_pen_combined, objs_true_combined, constraints_combined, heuristics_combined, designs_combined] = obtain_combined_data_case(truss_prob, model_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, constrad_prob_read, sidenodenum, n_pop, n_runs)
    
    n_members_total = nchoosek(sidenodenum^2,2);    
    data_array = zeros(n_pop,10,n_runs);
    if constrad_prob_read
        designs_array = strings(n_pop,n_runs);
        for i = 1:n_runs
            [data_array(:,:,i), designs_array(:,i)] = read_csv_data(truss_prob, model_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, constrad_prob_read, n_members_total, i-1);
        end
    else
        designs_array = zeros(n_pop,n_members_total,n_runs);
        for i = 1:n_runs
            [data_array(:,:,i), designs_array(:,:,i)] = read_csv_data(truss_prob, model_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, constrad_prob_read, n_members_total, i-1);
        end
    end
    [pen_obj1_combined, pen_obj2_combined, true_obj1_combined, true_obj2_combined, constraints_combined, heuristics_combined, designs_combined] = create_combined_arrays(data_array, designs_array, constrad_prob_read, n_members_total, n_runs);
    objs_pen_combined = [pen_obj1_combined, pen_obj2_combined];
    objs_true_combined = [true_obj1_combined, true_obj2_combined];
    
end

function [data_array, design_array] = read_csv_data(prob_truss, model_choice, partcoll_bools, nodalprop_bools, orient_bools, inters_bools, constrad_read, n_total_members, run_num)
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
        filepath_moea = "Epsilon MOEA\\";
    end
    
    switch model_choice
        case "Fibre"
            filepath3 = "Fibre Model\\";
            if truss_problem
                if constrad_read
                    filename2 = strcat(filename2,"_prob2_fibre.csv");
                else 
                    filename2 = strcat(filename2,"_fibre_varrad.csv");
                end
            else
                disp("Fiber stiffness model not suitable for artery problem")
                exit
            end
        case "Truss"
            filepath3 = "Truss Model\\";
            if prob_truss
                if constrad_read
                    filename2 = strcat(filename2,"_prob2_truss.csv");
                else
                    filename2 = strcat(filename2,"_truss_varrad.csv");
                end
            else
                filename2 = strcat(filename2,"_artery_truss.csv");
            end
        case "Beam"
            filepath3 = "Beam Model\\";
            if prob_truss
                if constrad_read
                    filename2 = strcat(filename2,"_prob2_beam.csv");
                else
                    filename2 = strcat(filename2,"_beam_varrad.csv");
                end
            else
                filename2 = strcat(filename2,"_artery_beam.csv");
            end
    end
    
    %%%% read appropriate file 
    full_filepath = strcat(filepath,filepath_prob,filepath_constrad,filepath3,filepath2,filepath_moea,filename,num2str(run_num),filename2);
    
    if prob_truss
        n_data = 11;
    else
        n_data = 10;
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

function [pen_obj1_combined, pen_obj2_combined, true_obj1_combined, true_obj2_combined, constr_combined, heur_combined, designs_combined] = create_combined_arrays(data_array, designs_array, read_constrad, n_total_members, n_runs)
    n_total = 0;
    n_constr = size(data_array(:,:,1),2) - 4 - 3; % number of constraints changes based on problem, so subtract 4 (2 pen. and 2 true objs.) and 3 (heurs) from number of total columns
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
    constr_combined = zeros(n_total,n_constr);
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
        constr_combined(index:index+n_current-1,:) = data_array_nonans(:,5:5+n_constr-1);
        heur_combined(index:index+n_current-1,:) = data_array_nonans(:,5+n_constr:10);
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

function top_edge_nodes = get_top_edge_nodes(sidenum)
    reached_right_edge = false;
    top_edge_nodes = zeros(sidenum,1);
    node = sidenum;
    count = 0;
    while(~reached_right_edge)
        if (node > (sidenum^2))
            reached_right_edge = true;
        else
            top_edge_nodes(count+1) = node;
            node = node + sidenum;
            count = count + 1;
        end
    end
end

function member_index = find_member_index(CA_full, member)
    member_index = 0;
    for i=1:size(CA_full,1)
        if (CA_full(i,1) == member(1))
            if (CA_full(i,2) == member(2))
                member_index = i;
            end
        end
    end
end

function norm_array = normalize_array(array)
    norm_array = zeros(size(array,1),size(array,2));
    for i = 1:size(array,2)
        array_col = array(:,i);
        max_array_col = max(array_col);
        min_array_col = min(array_col);
        norm_array_col = (array_col - min_array_col)/(max_array_col - min_array_col);
        norm_array(:,i) = norm_array_col;
    end
end

function objs_pf = compute_pareto_front(obj1_run, obj2_run)
    objs = [obj1_run, obj2_run];
    pf_bool = paretofront(objs);
    objs_pf = objs(pf_bool==1,:);
end

function min_dist_pf = compute_min_pf_dist(x_vals_norm, pf_objs_norm)
    dist_vals = zeros(size(pf_objs_norm,1),1);
    for i = 1:size(pf_objs_norm,1)
        dist_vals(i) = sqrt((x_vals_norm(1) - pf_objs_norm(i,1))^2 + (x_vals_norm(2) - pf_objs_norm(i,2))^2);
    end
    min_dist_pf = min(dist_vals);
end

