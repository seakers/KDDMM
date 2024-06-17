%% Heuristic ease of satisfaction tests and correlation tests (2D 3x3 constant radii)
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

%minimum_support = 0.005;

%% Cases to consider for GA data
constrad_read = true;
% Case 1 - Epsilon MOEA
repeat3x3_case1 = true;
case1_partcoll_bools = [false, false, false, false];
case1_nodalprop_bools = [false, false, false, false];
case1_orient_bools = [false, false, false, false];
case1_feas_bools = [false, false];

% Case 2 - AOS - Orientation 
%repeat3x3_case2 = true;
%case2_partcoll_bools = [false, false, false, false];
%case2_nodalprop_bools = [false, false, false, false];
%case2_orient_bools = [false, true, false, false];
%case2_feas_bools = [true, true];

%% Generate random architectures, objectives, constraints and heuristics
n_des = 100; % number of architectures to generate
n_runs = 10; % number of runs to generate "n_des" architectures

n_total_members = nchoosek(sidenum^2,2);
n_repeated_members = 2*nchoosek(sidenum,2);
n_variables = n_total_members - n_repeated_members;

add_ga_data = true;

bool_des_all = zeros(36,n_des,n_runs);
bool_des_map = containers.Map;

%orient_avg_all = zeros(n_des,n_runs);
orient_avg_norm_all = struct;

min_dist_pen_pf_all = struct;
min_dist_true_pf_all = struct;

c22_all = struct;
%c11_all = struct;
vol_frac_all = struct;
c22_true_all = struct;
vol_frac_true_all = struct;
c22_pen_all = struct;
vol_frac_pen_all = struct;

feas_all = struct;
stiff_rat_all = struct;
conn_all = struct;

coll_all = struct;
nod_all = struct;
%orient_all = zeros(n_des,n_runs);
orient_norm_all = struct;
inters_all = struct;

support_full_feas_runs = zeros(n_runs,1);
support_full_conn_runs = zeros(n_runs,1);
support_full_stiffrat_runs = zeros(n_runs,1);
support_full_coll_runs = zeros(n_runs,1);
support_full_nod_runs = zeros(n_runs,1);
support_full_orient_runs = zeros(n_runs,1);
support_full_inters_runs = zeros(n_runs,1);

n_pen_pf_runs = zeros(n_runs,1);
%n_true_pf_runs = zeros(n_runs,1);

support_pen_pf_runs = zeros(n_runs,1);
%support_true_pf_runs = zeros(n_runs,1);

if add_ga_data
    fibre_model_used = false;
    pop_size = 100;
    num_runs = 3;
    [f_pen_combined_case1, f_true_combined_case1, constr_combined_case1, heur_combined_case1, ~] = obtain_combined_data_case(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_feas_bools, constrad_read, repeat3x3_case1, sidenum, pop_size, num_runs);
    %[f_pen_combined_case2, f_true_combined_case2, constr_combined_case2, heur_combined_case2, ~] = obtain_combined_data_case(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_feas_bools, constrad_read, repeat3x3_case2, sidenum, pop_size, num_runs);
    
    %f_pen_combined_allcases = cat(1,f_pen_combined_case1,f_pen_combined_case2);
    %f_true_combined_allcases = cat(1,f_true_combined_case1,f_true_combined_case2);
    %constr_combined_allcases = cat(1,constr_combined_case1,constr_combined_case2);
    %heur_combined_allcases = cat(1,heur_combined_case1,heur_combined_case2);
    
    f_pen_combined_allcases = f_pen_combined_case1;
    f_true_combined_allcases = f_true_combined_case1;
    constr_combined_allcases = constr_combined_case1;
    heur_combined_allcases = heur_combined_case1;
    
    % constr = [feas,conn,stiffrat], heur = [partcoll,nodalprop,orient]
end
    
for i = 1:n_runs
    des_count = 1;
    %c11_run = zeros(n_des,1);
    c22_run = zeros(n_des,1);
    vol_frac_run = zeros(n_des,1);
    c22_pen_run = zeros(n_des,1);
    vol_frac_pen_run = zeros(n_des,1);
    feas_run = zeros(n_des,1);
    stiff_rat_run = zeros(n_des,1);
    conn_run = zeros(n_des,1);
    coll_run = zeros(n_des,1);
    nod_run = zeros(n_des,1);
    %orient_run = zeros(n_des,1);
    orient_norm_run = zeros(n_des,1);
    %inters_run = zeros(n_des,1);
    %orient_avg_run = zeros(n_des,1);
    %orient_avg_norm_run = zeros(n_des,1);
    bool_des_run = zeros(36,n_des);
    while des_count<(n_des+1)
        bool_des = randi([0,1],n_variables,1);
        %complete_bool_des = [bool_des(1:17); bool_des(3); bool_des(18:19); bool_des(6); bool_des(20:30); bool_des(22); bool_des(1); bool_des(2); bool_des(9)];
        complete_bool_des = get_complete_boolean_array(bool_des, CA_all, sidenum);
        CA_des = CA_all(complete_bool_des~=0,:);
        rvar_des = r.*ones(1,size(CA_des,1));
        [C_des, volfrac_des] = trussMetaCalc_NxN_rVar_AVar(nucFac,sel,rvar_des,E,CA_des);
        if (any(isnan(C_des),'all'))
            continue
        else
            if (map_contains_design(bool_des_map,bool_des))
                continue
            else
                bool_des_run(:,des_count) = complete_bool_des;
                bool_des_map(strcat('run',num2str(i),'des',num2str(des_count))) = bool_des_run;
                feas_run(des_count) = feasibility_checker_nonbinary_V2(NC,CA_des);
                %inters_run(des_count) = feas_run(des_count);
                %c11_run(des_count) = C_des(1,1);
                c22_run(des_count) = C_des(2,2);
                stiff_rat_run(des_count) = abs((C_des(2,2)/C_des(1,1)) - c_ratio);
                conn_run(des_count) = connectivityConstraint_NPBC_2D_V2(NC,CA_des,biasFactor);
                vol_frac_run(des_count) = volfrac_des;
                coll_run(des_count) = partCollapseHeuristic_2D(sidenum,CA_des,NC,sel,collapsibilityBiasFac);
                nod_run(des_count) = connectivityHeuristic_2D(sidenum,NC,CA_des,sel,biasFactor);
                %[orient_des, ~] = orientationHeuristic(NC,CA_des,c_ratio);
                [orient_norm_des, ~] = orientationHeuristicNorm(NC,CA_des,sel,c_ratio);
                %orient_run(des_count) = orient_des;
                orient_norm_run(des_count) = orient_norm_des;
                %orient_avg_run(des_count) = avg_angle_des;
                %orient_avg_norm_run(des_count) = avg_angle_norm_des;
                des_count = des_count + 1;
            end
        end
    end
    
    if add_ga_data
        c22_total = cat(1,c22_run,f_true_combined_allcases(:,1));
        volfrac_total = cat(1,vol_frac_run,f_true_combined_allcases(:,2));
        
        feas_total = cat(1,feas_run,constr_combined_allcases(:,1));
        conn_total = cat(1,conn_run,constr_combined_allcases(:,2));
        stiffrat_total = cat(1,stiff_rat_run,constr_combined_allcases(:,3));
        
        coll_total = cat(1,coll_run,heur_combined_allcases(:,1));
        nod_total = cat(1,nod_run,heur_combined_allcases(:,2));
        orient_total = cat(1,orient_norm_run,heur_combined_allcases(:,3));
        inters_total = feas_total;
    
    else
        c22_total = c22_run;
        volfrac_total = vol_frac_run;
        
        feas_total = feas_run;
        conn_total = conn_run;
        stiffrat_total = stiff_rat_run;
        
        coll_total = coll_run;
        nod_total = nod_run;
        orient_total = orient_norm_run;
        inters_total = feas_run;
    end
        
    current_field = strcat('trial',num2str(i));
    %c11_all(:,i) = c11_run;
    
    feas_all.(current_field) = feas_total;
    stiff_rat_all.(current_field) = stiffrat_total;
    conn_all.(current_field) = conn_total;
    
    stiffrat_constr_total = stiffrat_total./10;
    feas_pen_total = log10(abs(feas_total))./16;
    conn_pen_total = log10(abs(conn_total))./16;
    
    c22_true_all.(current_field) = c22_total;
    vol_frac_true_all.(current_field) = volfrac_total;
    
    c22_pen_total = -c22_total./E + (10.*stiffrat_constr_total) - 10.*(feas_pen_total + conn_pen_total)./2;  
    volfrac_pen_total = volfrac_total/0.96 + (10.*stiffrat_constr_total) - 10.*(feas_pen_total + conn_pen_total)./2;  
    
    pen_objs_pareto = compute_pareto_front(c22_pen_total,volfrac_pen_total);
    true_objs_pareto = compute_pareto_front(-c22_total,volfrac_total);
    true_objs_pareto_correct = [-true_objs_pareto(:,1),true_objs_pareto(:,2)];
    
    %%%% Normalizing objectives and pfs wrt max and min from objectives
    c22_max = max(c22_total);
    c22_min = min(c22_total);
    c22_pen_max = max(c22_pen_total);
    c22_pen_min = min(c22_pen_total);
    
    volfrac_max = max(volfrac_total);
    volfrac_min = min(volfrac_total);
    volfrac_pen_max = max(volfrac_pen_total);
    volfrac_pen_min = min(volfrac_pen_total);
    
    c22_pen_norm_total = (c22_pen_total - c22_pen_min)/(c22_pen_max - c22_pen_min);
    volfrac_pen_norm_total = (volfrac_pen_total - volfrac_pen_min)/(volfrac_pen_max - volfrac_pen_min);
    
    c22_norm_total = (c22_total - c22_min)/(c22_max - c22_min);
    volfrac_norm_total = (volfrac_total - volfrac_min)/(volfrac_max - volfrac_min);
    
    c22_pen_all.(current_field) = c22_pen_norm_total;
    vol_frac_pen_all.(current_field) = volfrac_pen_norm_total;
    
    c22_all.(current_field) = c22_norm_total;
    vol_frac_all.(current_field) = volfrac_norm_total;
    
    pen_objs_norm_pareto = [(pen_objs_pareto(:,1) - c22_pen_min)/(c22_pen_max - c22_pen_min), (pen_objs_pareto(:,2) - volfrac_pen_min)/(volfrac_pen_max - volfrac_pen_min)]; 
    true_objs_norm_pareto = [(true_objs_pareto_correct(:,1) - c22_min)/(c22_max - c22_min), (true_objs_pareto_correct(:,2) - volfrac_min)/(volfrac_max - volfrac_min)]; 
    
    min_dist_pen_pf_total = zeros(size(c22_pen_total,1),1);
    min_dist_true_pf_total = zeros(size(c22_pen_total,1),1);
    for k = 1:size(c22_pen_total,1)
        min_dist_pen_pf_total(k,1) = compute_min_pf_dist([c22_pen_norm_total(k,1),volfrac_pen_norm_total(k,1)],pen_objs_norm_pareto);
        min_dist_true_pf_total(k,1) = compute_min_pf_dist([c22_norm_total(k,1),volfrac_norm_total(k,1)],true_objs_norm_pareto);
    end
    
    bool_des_all(:,:,i) = bool_des_run;
    
    min_dist_pen_pf_all.(current_field) = min_dist_pen_pf_total;
    min_dist_true_pf_all.(current_field) = min_dist_true_pf_total;
    
    coll_all.(current_field) = coll_total;
    nod_all.(current_field) = nod_total;
    %orient_all(:,i) = orient_run;
    orient_norm_all.(current_field) = orient_total;
    inters_all.(current_field) = inters_total;
    %orient_avg_all(:,i) = orient_avg_run;
    %orient_avg_norm_all(:,i) = orient_avg_norm_run;
    
    support_full_feas_runs(i,1) = length(feas_total(feas_total==1))/size(c22_total,1);
    support_full_conn_runs(i,1) = length(conn_total(conn_total==1))/size(c22_total,1);
    support_full_stiffrat_runs(i,1) = length(stiffrat_total(stiffrat_total==0))/size(c22_total,1);
    support_full_coll_runs(i,1) = length(coll_total(coll_total==1))/size(c22_total,1);
    support_full_nod_runs(i,1) = length(nod_total(nod_total==1))/size(c22_total,1);
    support_full_orient_runs(i,1) = length(orient_total(orient_total==1))/size(c22_total,1);
    support_full_inters_runs(i,1) = length(inters_total(inters_total==1))/size(c22_total,1);
    
    n_pen_pf_runs(i,1) = size(pen_objs_pareto,1);
    
    support_pen_pf_runs(i,1) = size(pen_objs_pareto,1)/size(c22_total,1);

end
support_tablestats = [mean(support_pen_pf_runs), std(support_pen_pf_runs);
    mean(support_full_feas_runs), std(support_full_feas_runs);
    mean(support_full_conn_runs), std(support_full_conn_runs);
    mean(support_full_stiffrat_runs), std(support_full_stiffrat_runs);
    mean(support_full_coll_runs), std(support_full_coll_runs);
    mean(support_full_nod_runs), std(support_full_nod_runs);
    mean(support_full_orient_runs), std(support_full_orient_runs);
    mean(support_full_inters_runs), std(support_full_inters_runs)];

%% Heuristic violation plots
data_based_threshold = true;

n_des_total = 0;
for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    c22_total = c22_all.(current_field);
    n_des_total = n_des_total + size(c22_total,1);
end
c22_array = zeros(n_des_total,1);
volfrac_array = zeros(n_des_total,1);

c22_true_array = zeros(n_des_total,1);
volfrac_true_array = zeros(n_des_total,1);

feas_array = zeros(n_des_total,1);
conn_array = zeros(n_des_total,1);
stiffrat_array = zeros(n_des_total,1);

coll_array = zeros(n_des_total,1);
nod_array = zeros(n_des_total,1);
orient_array = zeros(n_des_total,1);
inters_array = zeros(n_des_total,1);

min_dist_pen_pf_array = zeros(n_des_total,1);
min_dist_true_pf_array = zeros(n_des_total,1);

index = 1;
for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    c22_total = c22_all.(current_field);
    volfrac_total = vol_frac_all.(current_field);
    
    c22_true_total = c22_true_all.(current_field);
    volfrac_true_total = vol_frac_true_all.(current_field);
    
    feas_total = feas_all.(current_field);
    conn_total = conn_all.(current_field);
    stiffrat_total = stiff_rat_all.(current_field);
    
    coll_total = coll_all.(current_field);
    nod_total = nod_all.(current_field);
    orient_total = orient_norm_all.(current_field);
    inters_total = inters_all.(current_field);
    
    min_dist_pen_pf_total = min_dist_pen_pf_all.(current_field);
    min_dist_true_pf_total = min_dist_true_pf_all.(current_field);
    
    n_des_run = size(c22_total,1);
    
    c22_array(index:index+n_des_run-1,1) = c22_total;
    volfrac_array(index:index+n_des_run-1,1) = volfrac_total;
    
    c22_true_array(index:index+n_des_run-1,1) = c22_true_total;
    volfrac_true_array(index:index+n_des_run-1,1) = volfrac_true_total;
    
    feas_array(index:index+n_des_run-1,1) = feas_total;
    conn_array(index:index+n_des_run-1,1) = conn_total;
    stiffrat_array(index:index+n_des_run-1,1) = stiffrat_total;
    
    coll_array(index:index+n_des_run-1,1) = coll_total;
    nod_array(index:index+n_des_run-1,1) = nod_total;
    orient_array(index:index+n_des_run-1,1) = orient_total;
    inters_array(index:index+n_des_run-1,1) = inters_total;
    
    min_dist_pen_pf_array(index:index+n_des_run-1,1) = min_dist_pen_pf_total;
    min_dist_true_pf_array(index:index+n_des_run-1,1) = min_dist_true_pf_total;

    index = index + n_des_run;
end

if data_based_threshold
    c22_thresh_val = prctile(c22_array,75);
    volfrac_thresh_val = prctile(volfrac_array,25);
    
    feas_thresh_val = prctile(feas_array,75);
    conn_thresh_val = prctile(conn_array,75);
    stiffrat_thresh_val = prctile(stiffrat_array,25);
    
    coll_thresh_val = prctile(coll_array,75);
    nod_thresh_val = prctile(nod_array,75);
    orient_thresh_val = prctile(orient_array,75);
    inters_thresh_val = prctile(inters_array,75);
    
    mindist_pfpen_thresh_val = prctile(min_dist_pen_pf_array,60);
    mindist_pftrue_thresh_val = prctile(min_dist_true_pf_array,70);
    
else
    c22_thresh_val = 0.7;
    volfrac_thresh_val = 0.3;
    
    feas_thresh_val = 0.5;
    conn_thresh_val = 0.9;
    stiffrat_thresh_val = 0.1;
    
    coll_thresh_val = 0.95;
    nod_thresh_val = 0.85;
    orient_thresh_val = 0.9;
    inters_thresh_val = 0.5;
    
    mindist_pfpen_thresh_val = 0.05;
    mindist_pftrue_thresh_val = 0.2;
end

stiffrat_constr_array = abs(stiffrat_array)/10;
feas_pen_array = log10(abs(feas_array))/16;
conn_pen_array = log10(abs(conn_array))/16;

c22_pen_array = -c22_true_array./E + (10.*stiffrat_constr_array) - 10.*(feas_pen_array + conn_pen_array)./2;  
volfrac_pen_array = volfrac_true_array/0.96 + (10.*stiffrat_constr_array) - 10.*(feas_pen_array + conn_pen_array)./2;  

score_tabelstats = [mean(feas_array), std(feas_array);
    mean(conn_array), std(conn_array);
    mean(stiffrat_array), std(stiffrat_array);
    mean(coll_array), std(coll_array);
    mean(nod_array), std(nod_array);
    mean(orient_array), std(orient_array);
    mean(inters_array), std(inters_array)];

figure 
scatter(c22_pen_array,volfrac_pen_array,[],1 - coll_array,'filled')
xlabel('Penalized $C_{22}$','Interpreter','Latex','FontSize',16)
ylabel('Penalized $v_f$','Interpreter','Latex','FontSize',16)
colorbar
title('Partial Collapsibility Violation','FontSize',16)

figure 
scatter(c22_pen_array,volfrac_pen_array,[],1 - nod_array,'filled')
xlabel('Penalized $C_{22}$','Interpreter','Latex','FontSize',16)
ylabel('Penalized $v_f$','Interpreter','Latex','FontSize',16)
colorbar;
title('Nodal Properties Violation','FontSize',16)

figure 
scatter(c22_pen_array,volfrac_pen_array,[],1 - orient_array,'filled')
xlabel('Penalized $C_{22}$','Interpreter','Latex','FontSize',16)
ylabel('Penalized $v_f$','Interpreter','Latex','FontSize',16)
colorbar;
title('Orientation Violation','FontSize',16)

figure 
scatter(c22_pen_array,volfrac_pen_array,[],1 - inters_array,'filled')
xlabel('Penalized $C_{22}$','Interpreter','Latex','FontSize',16)
ylabel('Penalized $v_f$','Interpreter','Latex','FontSize',16)
colorbar;
title('Intersection Violation','FontSize',16)

figure 
scatter(c22_true_array,volfrac_true_array,[],1 - coll_array,'filled')
xlabel('$C_{22}$','Interpreter','Latex')
ylabel('$v_f$','Interpreter','Latex')
colorbar
title('Partial Collapsibility Violation')

figure 
scatter(c22_true_array,volfrac_true_array,[],1 - nod_array,'filled')
xlabel('$C_{22}$','Interpreter','Latex')
ylabel('$v_f$','Interpreter','Latex')
colorbar;
title('Nodal Properties Violation')

figure 
scatter(c22_true_array,volfrac_true_array,[],1 - orient_array,'filled')
xlabel('$C_{22}$','Interpreter','Latex')
ylabel('$v_f$','Interpreter','Latex')
colorbar;
title('Orientation Violation')

figure 
scatter(c22_true_array,volfrac_true_array,[],1 - inters_array,'filled')
xlabel('$C_{22}$','Interpreter','Latex')
ylabel('$v_f$','Interpreter','Latex')
colorbar;
title('Intersection Violation')

%% Compute correlation coefficients of each heuristic with objectives and constraints
pearson_coll_truepfdist = zeros(n_runs,1);
pearson_coll_penpfdist = zeros(n_runs,1);
pearson_coll_feas = zeros(n_runs,1);
pearson_coll_stiffrat = zeros(n_runs,1);
pearson_coll_conn = zeros(n_runs,1);
% pval_coll_truepfdist = zeros(n_runs,1);
% pval_coll_feas = zeros(n_runs,1);
% pval_coll_stiffrat = zeros(n_runs,1);
% pval_coll_conn = zeros(n_runs,1);

pearson_nod_truepfdist = zeros(n_runs,1);
pearson_nod_penpfdist = zeros(n_runs,1);
pearson_nod_feas = zeros(n_runs,1);
pearson_nod_stiffrat = zeros(n_runs,1);
pearson_nod_conn = zeros(n_runs,1);
% pval_nod_truepfdist = zeros(n_runs,1);
% pval_nod_feas = zeros(n_runs,1);
% pval_nod_stiffrat = zeros(n_runs,1);
% pval_nod_conn = zeros(n_runs,1);

pearson_orient_truepfdist = zeros(n_runs,1);
pearson_orient_penpfdist = zeros(n_runs,1);
pearson_orient_feas = zeros(n_runs,1);
pearson_orient_stiffrat = zeros(n_runs,1);
pearson_orient_conn = zeros(n_runs,1);
% pval_orient_truepfdist = zeros(n_runs,1);
% pval_orient_feas = zeros(n_runs,1);
% pval_orient_stiffrat = zeros(n_runs,1);
% pval_orient_conn = zeros(n_runs,1);

pearson_inters_truepfdist = zeros(n_runs,1);
pearson_inters_penpfdist = zeros(n_runs,1);
pearson_inters_feas = zeros(n_runs,1);
pearson_inters_stiffrat = zeros(n_runs,1);
pearson_inters_conn = zeros(n_runs,1);
% pval_inters_truepfdist = zeros(n_runs,1);
% pval_inters_feas = zeros(n_runs,1);
% pval_inters_stiffrat = zeros(n_runs,1);
% pval_inters_conn = zeros(n_runs,1);

spearman_coll_truepfdist = zeros(n_runs,1);
spearman_coll_penpfdist = zeros(n_runs,1);
spearman_coll_feas = zeros(n_runs,1);
spearman_coll_stiffrat = zeros(n_runs,1);
spearman_coll_conn = zeros(n_runs,1);
% pval_coll_truepfdist = zeros(n_runs,1);
% pval_coll_feas = zeros(n_runs,1);
% pval_coll_stiffrat = zeros(n_runs,1);
% pval_coll_conn = zeros(n_runs,1);

spearman_nod_truepfdist = zeros(n_runs,1);
spearman_nod_penpfdist = zeros(n_runs,1);
spearman_nod_feas = zeros(n_runs,1);
spearman_nod_stiffrat = zeros(n_runs,1);
spearman_nod_conn = zeros(n_runs,1);
% pval_nod_truepfdist = zeros(n_runs,1);
% pval_nod_feas = zeros(n_runs,1);
% pval_nod_stiffrat = zeros(n_runs,1);
% pval_nod_conn = zeros(n_runs,1);

spearman_orient_truepfdist = zeros(n_runs,1);
spearman_orient_penpfdist = zeros(n_runs,1);
spearman_orient_feas = zeros(n_runs,1);
spearman_orient_stiffrat = zeros(n_runs,1);
spearman_orient_conn = zeros(n_runs,1);
% pval_orient_truepfdist = zeros(n_runs,1);
% pval_orient_feas = zeros(n_runs,1);
% pval_orient_stiffrat = zeros(n_runs,1);
% pval_orient_conn = zeros(n_runs,1);

spearman_inters_truepfdist = zeros(n_runs,1);
spearman_inters_penpfdist = zeros(n_runs,1);
spearman_inters_feas = zeros(n_runs,1);
spearman_inters_stiffrat = zeros(n_runs,1);
spearman_inters_conn = zeros(n_runs,1);
% pval_inters_truepfdist = zeros(n_runs,1);
% pval_inters_feas = zeros(n_runs,1);
% pval_inters_stiffrat = zeros(n_runs,1);
% pval_inters_conn = zeros(n_runs,1);

for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    min_dist_true_pf_total = min_dist_true_pf_all.(current_field);
    min_dist_pen_pf_total = min_dist_pen_pf_all.(current_field);
    
    coll_total = coll_all.(current_field);
    nod_total = nod_all.(current_field);
    orient_total = orient_norm_all.(current_field);
    inters_total = inters_all.(current_field);
    
    feas_total = feas_all.(current_field);
    conn_total = conn_all.(current_field);
    stiffrat_total = stiff_rat_all.(current_field);
    
    [pearson_coll_truepfdist(i),~] = corr(coll_total,min_dist_true_pf_total,'Type','Pearson','Rows','complete');
    [pearson_coll_penpfdist(i),~] = corr(coll_total,min_dist_pen_pf_total,'Type','Pearson','Rows','complete');
    [pearson_coll_feas(i),~] = corr(coll_total,feas_total,'Type','Pearson','Rows','complete');
    [pearson_coll_stiffrat(i),~] = corr((coll_total),stiffrat_total,'Type','Pearson','Rows','complete');
    [pearson_coll_conn(i),~] = corr(coll_total,conn_total,'Type','Pearson','Rows','complete');
    
    [pearson_nod_truepfdist(i),~] = corr(nod_total,min_dist_true_pf_total,'Type','Pearson','Rows','complete');
    [pearson_nod_penpfdist(i),~] = corr(nod_total,min_dist_pen_pf_total,'Type','Pearson','Rows','complete');
    [pearson_nod_feas(i),~] = corr(nod_total,feas_total,'Type','Pearson','Rows','complete');
    [pearson_nod_stiffrat(i),~] = corr(nod_total,stiffrat_total,'Type','Pearson','Rows','complete');
    [pearson_nod_conn(i),~] = corr(nod_total,conn_total,'Type','Pearson','Rows','complete');
    
    [pearson_orient_truepfdist(i),~] = corr(orient_total,min_dist_true_pf_total,'Type','Pearson','Rows','complete');
    [pearson_orient_penpfdist(i),~] = corr(orient_total,min_dist_pen_pf_total,'Type','Pearson','Rows','complete');
    [pearson_orient_feas(i),~] = corr(orient_total,feas_total,'Type','Pearson','Rows','complete');
    [pearson_orient_stiffrat(i),~] = corr(orient_total,stiffrat_total,'Type','Pearson','Rows','complete');
    [pearson_orient_conn(i),~] = corr(orient_total,conn_total,'Type','Pearson','Rows','complete');
    
    [pearson_inters_truepfdist(i),~] = corr(inters_total,min_dist_true_pf_total,'Type','Pearson','Rows','complete');
    [pearson_inters_penpfdist(i),~] = corr(inters_total,min_dist_pen_pf_total,'Type','Pearson','Rows','complete');
    [pearson_inters_feas(i),~] = corr(inters_total,feas_total,'Type','Pearson','Rows','complete');
    [pearson_inters_stiffrat(i),~] = corr(inters_total,stiffrat_total,'Type','Pearson','Rows','complete');
    [pearson_inters_conn(i),~] = corr(inters_total,conn_total,'Type','Pearson','Rows','complete');
    
    [spearman_coll_truepfdist(i),~] = corr(coll_total,min_dist_true_pf_total,'Type','Spearman','Rows','complete');
    [spearman_coll_penpfdist(i),~] = corr(coll_total,min_dist_pen_pf_total,'Type','Spearman','Rows','complete');
    [spearman_coll_feas(i),~] = corr(coll_total,feas_total,'Type','Spearman','Rows','complete');
    [spearman_coll_stiffrat(i),~] = corr((coll_total),stiffrat_total,'Type','Spearman','Rows','complete');
    [spearman_coll_conn(i),~] = corr(coll_total,conn_total,'Type','Spearman','Rows','complete');
    
    [spearman_nod_truepfdist(i),~] = corr(nod_total,min_dist_true_pf_total,'Type','Spearman','Rows','complete');
    [spearman_nod_penpfdist(i),~] = corr(nod_total,min_dist_pen_pf_total,'Type','Spearman','Rows','complete');
    [spearman_nod_feas(i),~] = corr(nod_total,feas_total,'Type','Spearman','Rows','complete');
    [spearman_nod_stiffrat(i),~] = corr(nod_total,stiffrat_total,'Type','Spearman','Rows','complete');
    [spearman_nod_conn(i),~] = corr(nod_total,conn_total,'Type','Spearman','Rows','complete');
    
    [spearman_orient_truepfdist(i),~] = corr(orient_total,min_dist_true_pf_total,'Type','Spearman','Rows','complete');
    [spearman_orient_penpfdist(i),~] = corr(orient_total,min_dist_pen_pf_total,'Type','Spearman','Rows','complete');
    [spearman_orient_feas(i),~] = corr(orient_total,feas_total,'Type','Spearman','Rows','complete');
    [spearman_orient_stiffrat(i),~] = corr(orient_total,stiffrat_total,'Type','Spearman','Rows','complete');
    [spearman_orient_conn(i),~] = corr(orient_total,conn_total,'Type','Spearman','Rows','complete');
    
    [spearman_inters_truepfdist(i),~] = corr(inters_total,min_dist_true_pf_total,'Type','Spearman','Rows','complete');
    [spearman_inters_penpfdist(i),~] = corr(inters_total,min_dist_pen_pf_total,'Type','Spearman','Rows','complete');
    [spearman_inters_feas(i),~] = corr(inters_total,feas_total,'Type','Spearman','Rows','complete');
    [spearman_inters_stiffrat(i),~] = corr(inters_total,stiffrat_total,'Type','Spearman','Rows','complete');
    [spearman_inters_conn(i),~] = corr(inters_total,conn_total,'Type','Spearman','Rows','complete');
    
end

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
   
correlation_penpf_tablestats = [mean(pearson_coll_penpfdist), std(pearson_coll_penpfdist), mean(spearman_coll_penpfdist), std(spearman_coll_penpfdist);
    mean(pearson_nod_penpfdist), std(pearson_nod_penpfdist), mean(spearman_nod_penpfdist), std(spearman_nod_penpfdist);
    mean(pearson_orient_penpfdist), std(pearson_orient_penpfdist), mean(spearman_orient_penpfdist), std(spearman_orient_penpfdist);
    mean(pearson_inters_penpfdist), std(pearson_inters_penpfdist), mean(spearman_inters_penpfdist), std(spearman_inters_penpfdist)];
    
%% Thresholding heuristics, objectives and constraints into high and low
coll_all_thresh = struct;
nod_all_thresh = struct;
orient_all_thresh = struct;
inters_all_thresh = struct;

c22_all_thresh = struct;
volfrac_all_thresh = struct;

feas_all_thresh = struct;
conn_all_thresh = struct;
stiffrat_all_thresh = struct;

min_dist_pen_pf_all_thresh = struct;
min_dist_true_pf_all_thresh = struct;

for i = 1:n_runs
    current_field = strcat('trial',num2str(i));
    coll_total = coll_all.(current_field);
    nod_total = nod_all.(current_field);
    orient_total = orient_norm_all.(current_field);
    inters_total = inters_all.(current_field);
    
    c22_total = c22_all.(current_field);
    volfrac_total = vol_frac_all.(current_field);
    
    feas_total = feas_all.(current_field);
    stiffrat_total = stiff_rat_all.(current_field);
    conn_total = conn_all.(current_field);
    
    min_dist_pen_pf_total = min_dist_pen_pf_all.(current_field);
    min_dist_true_pf_total = min_dist_true_pf_all.(current_field);
    
    coll_run_thresh = double(coll_total >= coll_thresh_val);
    nod_run_thresh = double(nod_total >= nod_thresh_val);
    orient_run_thresh = double(orient_total >= orient_thresh_val);
    inters_run_thresh = double(inters_total >= inters_thresh_val);
    
    c22_run_thresh = double(c22_total >= c22_thresh_val);
    volfrac_run_thresh = double(volfrac_total >= volfrac_thresh_val);
    
    feas_run_thresh = double(feas_total >= feas_thresh_val);
    conn_run_thresh = double(conn_total >= conn_thresh_val);
    stiffrat_run_thresh = double(stiffrat_total >= stiffrat_thresh_val);
    
    min_dist_pen_pf_thresh = double(min_dist_pen_pf_total >= mindist_pfpen_thresh_val);
    min_dist_true_pf_thresh = double(min_dist_true_pf_total >= mindist_pftrue_thresh_val);
    
    coll_all_thresh.(current_field) = coll_run_thresh;
    nod_all_thresh.(current_field) = nod_run_thresh;
    orient_all_thresh.(current_field) = orient_run_thresh;
    inters_all_thresh.(current_field) = inters_run_thresh;
    
    c22_all_thresh.(current_field) = c22_run_thresh;
    volfrac_all_thresh.(current_field) = volfrac_run_thresh;
    
    feas_all_thresh.(current_field) = feas_run_thresh;
    conn_all_thresh.(current_field) = conn_run_thresh;
    stiffrat_all_thresh.(current_field) = stiffrat_run_thresh;
    
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
    c22_thresh = c22_all_thresh.(current_field);
    volfrac_thresh = volfrac_all_thresh.(current_field);
    
    coll_thresh = coll_all_thresh.(current_field);
    nod_thresh = nod_all_thresh.(current_field);
    orient_thresh = orient_all_thresh.(current_field);
    inters_thresh = inters_all_thresh.(current_field);
    
    feas_thresh = feas_all_thresh.(current_field);
    conn_thresh = conn_all_thresh.(current_field);
    stiffrat_thresh = stiffrat_all_thresh.(current_field);
    
    min_dist_pen_pf_thresh = min_dist_pen_pf_all_thresh.(current_field);
    min_dist_true_pf_thresh = min_dist_true_pf_all_thresh.(current_field);
    
    coll_intrness_truepf_run = compute_heur_intrness_penpf(coll_thresh, min_dist_true_pf_thresh);
    coll_intrness_truepf.(current_field) = coll_intrness_truepf_run;
    nod_intrness_truepf_run = compute_heur_intrness_penpf(nod_thresh, min_dist_true_pf_thresh);
    nod_intrness_truepf.(current_field) = nod_intrness_truepf_run;
    orient_intrness_truepf_run = compute_heur_intrness_penpf(orient_thresh, min_dist_true_pf_thresh);
    orient_intrness_truepf.(current_field) = orient_intrness_truepf_run;
    inters_intrness_truepf_run = compute_heur_intrness_penpf(inters_thresh, min_dist_true_pf_thresh);
    inters_intrness_truepf.(current_field) = inters_intrness_truepf_run;
    
    coll_intrness_run = compute_heur_intrness_array(coll_thresh, min_dist_true_pf_thresh, c22_thresh, volfrac_thresh, feas_thresh, conn_thresh, stiffrat_thresh);
    coll_intrness_truepf_constr.(current_field) = coll_intrness_run;
    nod_intrness_run = compute_heur_intrness_array(nod_thresh, min_dist_true_pf_thresh, c22_thresh, volfrac_thresh, feas_thresh, conn_thresh, stiffrat_thresh);
    nod_intrness_truepf_constr.(current_field) = nod_intrness_run;
    orient_intrness_run = compute_heur_intrness_array(orient_thresh, min_dist_true_pf_thresh, c22_thresh, volfrac_thresh, feas_thresh, conn_thresh, stiffrat_thresh);
    orient_intrness_truepf_constr.(current_field) = orient_intrness_run;
    inters_intrness_run = compute_heur_intrness_array(inters_thresh, min_dist_true_pf_thresh, c22_thresh, volfrac_thresh, feas_thresh, conn_thresh, stiffrat_thresh);
    inters_intrness_truepf_constr.(current_field) = inters_intrness_run;
    
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
low_heur_close2truepf_lowc22_runs = zeros(n_runs,6);
low_heur_close2truepf_highc22_runs = zeros(n_runs,6);
low_heur_close2truepf_lowvf_runs = zeros(n_runs,6);
low_heur_close2truepf_highvf_runs = zeros(n_runs,6);
low_heur_low_feas_runs = zeros(n_runs,6);
low_heur_high_feas_runs = zeros(n_runs,6);
low_heur_low_conn_runs = zeros(n_runs,6);
low_heur_high_conn_runs = zeros(n_runs,6);
low_heur_low_stiffrat_runs = zeros(n_runs,6);
low_heur_high_stiffrat_runs = zeros(n_runs,6);

high_heur_close2truepf_lowc22_runs = zeros(n_runs,6);
high_heur_close2truepf_highc22_runs = zeros(n_runs,6);
high_heur_close2truepf_lowvf_runs = zeros(n_runs,6);
high_heur_close2truepf_highvf_runs = zeros(n_runs,6);
high_heur_low_feas_runs = zeros(n_runs,6);
high_heur_high_feas_runs = zeros(n_runs,6);
high_heur_low_conn_runs = zeros(n_runs,6);
high_heur_high_conn_runs = zeros(n_runs,6);
high_heur_low_stiffrat_runs = zeros(n_runs,6);
high_heur_high_stiffrat_runs = zeros(n_runs,6);

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
    
    low_heur_close2truepf_lowc22_runs(i,:) = heur_intrness_run(1,:,1);
    low_heur_close2truepf_highc22_runs(i,:) = heur_intrness_run(2,:,1);
    low_heur_close2truepf_lowvf_runs(i,:) = heur_intrness_run(3,:,1);
    low_heur_close2truepf_highvf_runs(i,:) = heur_intrness_run(4,:,1);
    low_heur_low_feas_runs(i,:) = heur_intrness_run(5,:,1);
    low_heur_high_feas_runs(i,:) = heur_intrness_run(6,:,1);
    low_heur_low_conn_runs(i,:) = heur_intrness_run(7,:,1);
    low_heur_high_conn_runs(i,:) = heur_intrness_run(8,:,1);
    low_heur_low_stiffrat_runs(i,:) = heur_intrness_run(9,:,1);
    low_heur_high_stiffrat_runs(i,:) = heur_intrness_run(10,:,1);
    
    high_heur_close2truepf_lowc22_runs(i,:) = heur_intrness_run(1,:,2);
    high_heur_close2truepf_highc22_runs(i,:) = heur_intrness_run(2,:,2);
    high_heur_close2truepf_lowvf_runs(i,:) = heur_intrness_run(3,:,2);
    high_heur_close2truepf_highvf_runs(i,:) = heur_intrness_run(4,:,2);
    high_heur_low_feas_runs(i,:) = heur_intrness_run(5,:,2);
    high_heur_high_feas_runs(i,:) = heur_intrness_run(6,:,2);
    high_heur_low_conn_runs(i,:) = heur_intrness_run(7,:,2);
    high_heur_high_conn_runs(i,:) = heur_intrness_run(8,:,2);
    high_heur_low_stiffrat_runs(i,:) = heur_intrness_run(9,:,2);
    high_heur_high_stiffrat_runs(i,:) = heur_intrness_run(10,:,2);
    
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

low_heur_tablestats = [mean(low_heur_close2truepf_lowc22_runs(:,1)), std(low_heur_close2truepf_lowc22_runs(:,1)), mean(low_heur_close2truepf_lowc22_runs(:,2)), std(low_heur_close2truepf_lowc22_runs(:,2)), mean(low_heur_close2truepf_lowc22_runs(:,3)), std(low_heur_close2truepf_lowc22_runs(:,3)), mean(low_heur_close2truepf_lowc22_runs(:,4)), std(low_heur_close2truepf_lowc22_runs(:,4)), mean(low_heur_close2truepf_lowc22_runs(:,5)), std(low_heur_close2truepf_lowc22_runs(:,5)), mean(low_heur_close2truepf_lowc22_runs(:,6)), std(low_heur_close2truepf_lowc22_runs(:,6));
    mean(low_heur_close2truepf_highc22_runs(:,1)), std(low_heur_close2truepf_highc22_runs(:,1)), mean(low_heur_close2truepf_highc22_runs(:,2)), std(low_heur_close2truepf_highc22_runs(:,2)), mean(low_heur_close2truepf_highc22_runs(:,3)), std(low_heur_close2truepf_highc22_runs(:,3)), mean(low_heur_close2truepf_highc22_runs(:,4)), std(low_heur_close2truepf_highc22_runs(:,4)), mean(low_heur_close2truepf_highc22_runs(:,5)), std(low_heur_close2truepf_highc22_runs(:,5)), mean(low_heur_close2truepf_highc22_runs(:,6)), std(low_heur_close2truepf_highc22_runs(:,6));
    mean(low_heur_close2truepf_lowvf_runs(:,1)), std(low_heur_close2truepf_lowvf_runs(:,1)), mean(low_heur_close2truepf_lowvf_runs(:,2)), std(low_heur_close2truepf_lowvf_runs(:,2)), mean(low_heur_close2truepf_lowvf_runs(:,3)), std(low_heur_close2truepf_lowvf_runs(:,3)), mean(low_heur_close2truepf_lowvf_runs(:,4)), std(low_heur_close2truepf_lowvf_runs(:,4)), mean(low_heur_close2truepf_lowvf_runs(:,5)), std(low_heur_close2truepf_lowvf_runs(:,5)), mean(low_heur_close2truepf_lowvf_runs(:,6)), std(low_heur_close2truepf_lowvf_runs(:,6));
    mean(low_heur_close2truepf_highvf_runs(:,1)), std(low_heur_close2truepf_highvf_runs(:,1)), mean(low_heur_close2truepf_highvf_runs(:,2)), std(low_heur_close2truepf_highvf_runs(:,2)), mean(low_heur_close2truepf_highvf_runs(:,3)), std(low_heur_close2truepf_highvf_runs(:,3)), mean(low_heur_close2truepf_highvf_runs(:,4)), std(low_heur_close2truepf_highvf_runs(:,4)), mean(low_heur_close2truepf_highvf_runs(:,5)), std(low_heur_close2truepf_highvf_runs(:,5)), mean(low_heur_close2truepf_highvf_runs(:,6)), std(low_heur_close2truepf_highvf_runs(:,6));
    mean(low_heur_low_feas_runs(:,1)), std(low_heur_low_feas_runs(:,1)), mean(low_heur_low_feas_runs(:,2)), std(low_heur_low_feas_runs(:,2)), mean(low_heur_low_feas_runs(:,3)), std(low_heur_low_feas_runs(:,3)), mean(low_heur_low_feas_runs(:,4)), std(low_heur_low_feas_runs(:,4)), mean(low_heur_low_feas_runs(:,5)), std(low_heur_low_feas_runs(:,5)), mean(low_heur_low_feas_runs(:,6)), std(low_heur_low_feas_runs(:,6));
    mean(low_heur_high_feas_runs(:,1)), std(low_heur_high_feas_runs(:,1)), mean(low_heur_high_feas_runs(:,2)), std(low_heur_high_feas_runs(:,2)), mean(low_heur_high_feas_runs(:,3)), std(low_heur_high_feas_runs(:,3)), mean(low_heur_high_feas_runs(:,4)), std(low_heur_high_feas_runs(:,4)), mean(low_heur_high_feas_runs(:,5)), std(low_heur_high_feas_runs(:,5)), mean(low_heur_high_feas_runs(:,6)), std(low_heur_high_feas_runs(:,6));
    mean(low_heur_low_conn_runs(:,1)), std(low_heur_low_conn_runs(:,1)), mean(low_heur_low_conn_runs(:,2)), std(low_heur_low_conn_runs(:,2)), mean(low_heur_low_conn_runs(:,3)), std(low_heur_low_conn_runs(:,3)), mean(low_heur_low_conn_runs(:,4)), std(low_heur_low_conn_runs(:,4)), mean(low_heur_low_conn_runs(:,5)), std(low_heur_low_conn_runs(:,5)), mean(low_heur_low_conn_runs(:,6)), std(low_heur_low_conn_runs(:,6));
    mean(low_heur_high_conn_runs(:,1)), std(low_heur_high_conn_runs(:,1)), mean(low_heur_high_conn_runs(:,2)), std(low_heur_high_conn_runs(:,2)), mean(low_heur_high_conn_runs(:,3)), std(low_heur_high_conn_runs(:,3)), mean(low_heur_high_conn_runs(:,4)), std(low_heur_high_conn_runs(:,4)), mean(low_heur_high_conn_runs(:,5)), std(low_heur_high_conn_runs(:,5)), mean(low_heur_high_conn_runs(:,6)), std(low_heur_high_conn_runs(:,6));
    mean(low_heur_low_stiffrat_runs(:,1)), std(low_heur_low_stiffrat_runs(:,1)), mean(low_heur_low_stiffrat_runs(:,2)), std(low_heur_low_stiffrat_runs(:,2)), mean(low_heur_low_stiffrat_runs(:,3)), std(low_heur_low_stiffrat_runs(:,3)), mean(low_heur_low_stiffrat_runs(:,4)), std(low_heur_low_stiffrat_runs(:,4)), mean(low_heur_low_stiffrat_runs(:,5)), std(low_heur_low_stiffrat_runs(:,5)), mean(low_heur_low_stiffrat_runs(:,6)), std(low_heur_low_stiffrat_runs(:,6));
    mean(low_heur_high_stiffrat_runs(:,1)), std(low_heur_high_stiffrat_runs(:,1)), mean(low_heur_high_stiffrat_runs(:,2)), std(low_heur_high_stiffrat_runs(:,2)), mean(low_heur_high_stiffrat_runs(:,3)), std(low_heur_high_stiffrat_runs(:,3)), mean(low_heur_high_stiffrat_runs(:,4)), std(low_heur_high_stiffrat_runs(:,4)), mean(low_heur_high_stiffrat_runs(:,5)), std(low_heur_high_stiffrat_runs(:,5)), mean(low_heur_high_stiffrat_runs(:,6)), std(low_heur_high_stiffrat_runs(:,6))];

high_heur_tablestats = [mean(high_heur_close2truepf_lowc22_runs(:,1)), std(high_heur_close2truepf_lowc22_runs(:,1)), mean(high_heur_close2truepf_lowc22_runs(:,2)), std(high_heur_close2truepf_lowc22_runs(:,2)), mean(high_heur_close2truepf_lowc22_runs(:,3)), std(high_heur_close2truepf_lowc22_runs(:,3)), mean(high_heur_close2truepf_lowc22_runs(:,4)), std(high_heur_close2truepf_lowc22_runs(:,4)), mean(high_heur_close2truepf_lowc22_runs(:,5)), std(high_heur_close2truepf_lowc22_runs(:,5)), mean(high_heur_close2truepf_lowc22_runs(:,6)), std(high_heur_close2truepf_lowc22_runs(:,6));
    mean(high_heur_close2truepf_highc22_runs(:,1)), std(high_heur_close2truepf_highc22_runs(:,1)), mean(high_heur_close2truepf_highc22_runs(:,2)), std(high_heur_close2truepf_highc22_runs(:,2)), mean(high_heur_close2truepf_highc22_runs(:,3)), std(high_heur_close2truepf_highc22_runs(:,3)), mean(high_heur_close2truepf_highc22_runs(:,4)), std(high_heur_close2truepf_highc22_runs(:,4)), mean(high_heur_close2truepf_highc22_runs(:,5)), std(high_heur_close2truepf_highc22_runs(:,5)), mean(high_heur_close2truepf_highc22_runs(:,6)), std(high_heur_close2truepf_highc22_runs(:,6));
    mean(high_heur_close2truepf_lowvf_runs(:,1)), std(high_heur_close2truepf_lowvf_runs(:,1)), mean(high_heur_close2truepf_lowvf_runs(:,2)), std(high_heur_close2truepf_lowvf_runs(:,2)), mean(high_heur_close2truepf_lowvf_runs(:,3)), std(high_heur_close2truepf_lowvf_runs(:,3)), mean(high_heur_close2truepf_lowvf_runs(:,4)), std(high_heur_close2truepf_lowvf_runs(:,4)), mean(high_heur_close2truepf_lowvf_runs(:,5)), std(high_heur_close2truepf_lowvf_runs(:,5)), mean(high_heur_close2truepf_lowvf_runs(:,6)), std(high_heur_close2truepf_lowvf_runs(:,6));
    mean(high_heur_close2truepf_highvf_runs(:,1)), std(high_heur_close2truepf_highvf_runs(:,1)), mean(high_heur_close2truepf_highvf_runs(:,2)), std(high_heur_close2truepf_highvf_runs(:,2)), mean(high_heur_close2truepf_highvf_runs(:,3)), std(high_heur_close2truepf_highvf_runs(:,3)), mean(high_heur_close2truepf_highvf_runs(:,4)), std(high_heur_close2truepf_highvf_runs(:,4)), mean(high_heur_close2truepf_highvf_runs(:,5)), std(high_heur_close2truepf_highvf_runs(:,5)), mean(high_heur_close2truepf_highvf_runs(:,6)), std(high_heur_close2truepf_highvf_runs(:,6));
    mean(high_heur_low_feas_runs(:,1)), std(high_heur_low_feas_runs(:,1)), mean(high_heur_low_feas_runs(:,2)), std(high_heur_low_feas_runs(:,2)), mean(high_heur_low_feas_runs(:,3)), std(high_heur_low_feas_runs(:,3)), mean(high_heur_low_feas_runs(:,4)), std(high_heur_low_feas_runs(:,4)), mean(high_heur_low_feas_runs(:,5)), std(high_heur_low_feas_runs(:,5)), mean(high_heur_low_feas_runs(:,6)), std(high_heur_low_feas_runs(:,6));
    mean(high_heur_high_feas_runs(:,1)), std(high_heur_high_feas_runs(:,1)), mean(high_heur_high_feas_runs(:,2)), std(high_heur_high_feas_runs(:,2)), mean(high_heur_high_feas_runs(:,3)), std(high_heur_high_feas_runs(:,3)), mean(high_heur_high_feas_runs(:,4)), std(high_heur_high_feas_runs(:,4)), mean(high_heur_high_feas_runs(:,5)), std(high_heur_high_feas_runs(:,5)), mean(high_heur_high_feas_runs(:,6)), std(high_heur_high_feas_runs(:,6));
    mean(high_heur_low_conn_runs(:,1)), std(high_heur_low_conn_runs(:,1)), mean(high_heur_low_conn_runs(:,2)), std(high_heur_low_conn_runs(:,2)), mean(high_heur_low_conn_runs(:,3)), std(high_heur_low_conn_runs(:,3)), mean(high_heur_low_conn_runs(:,4)), std(high_heur_low_conn_runs(:,4)), mean(high_heur_low_conn_runs(:,5)), std(high_heur_low_conn_runs(:,5)), mean(high_heur_low_conn_runs(:,6)), std(high_heur_low_conn_runs(:,6));
    mean(high_heur_high_conn_runs(:,1)), std(high_heur_high_conn_runs(:,1)), mean(high_heur_high_conn_runs(:,2)), std(high_heur_high_conn_runs(:,2)), mean(high_heur_high_conn_runs(:,3)), std(high_heur_high_conn_runs(:,3)), mean(high_heur_high_conn_runs(:,4)), std(high_heur_high_conn_runs(:,4)), mean(high_heur_high_conn_runs(:,5)), std(high_heur_high_conn_runs(:,5)), mean(high_heur_high_conn_runs(:,6)), std(high_heur_high_conn_runs(:,6));
    mean(high_heur_low_stiffrat_runs(:,1)), std(high_heur_low_stiffrat_runs(:,1)), mean(high_heur_low_stiffrat_runs(:,2)), std(high_heur_low_stiffrat_runs(:,2)), mean(high_heur_low_stiffrat_runs(:,3)), std(high_heur_low_stiffrat_runs(:,3)), mean(high_heur_low_stiffrat_runs(:,4)), std(high_heur_low_stiffrat_runs(:,4)), mean(high_heur_low_stiffrat_runs(:,5)), std(high_heur_low_stiffrat_runs(:,5)), mean(high_heur_low_stiffrat_runs(:,6)), std(high_heur_low_stiffrat_runs(:,6));
    mean(high_heur_high_stiffrat_runs(:,1)), std(high_heur_high_stiffrat_runs(:,1)), mean(high_heur_high_stiffrat_runs(:,2)), std(high_heur_high_stiffrat_runs(:,2)), mean(high_heur_high_stiffrat_runs(:,3)), std(high_heur_high_stiffrat_runs(:,3)), mean(high_heur_high_stiffrat_runs(:,4)), std(high_heur_high_stiffrat_runs(:,4)), mean(high_heur_high_stiffrat_runs(:,5)), std(high_heur_high_stiffrat_runs(:,5)), mean(high_heur_high_stiffrat_runs(:,6)), std(high_heur_high_stiffrat_runs(:,6))];

heur_penpf_tablestats = [mean(low_heur_close2penpf(:,1)), std(low_heur_close2penpf(:,1)), mean(low_heur_close2penpf(:,2)), std(low_heur_close2penpf(:,2)), mean(low_heur_close2penpf(:,3)), std(low_heur_close2penpf(:,3)), mean(low_heur_close2penpf(:,4)), std(low_heur_close2penpf(:,4)), mean(low_heur_close2penpf(:,5)), std(low_heur_close2penpf(:,5)), mean(low_heur_close2penpf(:,6)), std(low_heur_close2penpf(:,6));
    mean(high_heur_close2penpf(:,1)), std(high_heur_close2penpf(:,1)), mean(high_heur_close2penpf(:,2)), std(high_heur_close2penpf(:,2)), mean(high_heur_close2penpf(:,3)), std(high_heur_close2penpf(:,3)), mean(high_heur_close2penpf(:,4)), std(high_heur_close2penpf(:,4)), mean(high_heur_close2penpf(:,5)), std(high_heur_close2penpf(:,5)), mean(high_heur_close2penpf(:,6)), std(high_heur_close2penpf(:,6))];

heur_truepf_tablestats = [mean(low_heur_close2truepf(:,1)), std(low_heur_close2truepf(:,1)), mean(low_heur_close2truepf(:,2)), std(low_heur_close2truepf(:,2)), mean(low_heur_close2truepf(:,3)), std(low_heur_close2truepf(:,3)), mean(low_heur_close2truepf(:,4)), std(low_heur_close2truepf(:,4)), mean(low_heur_close2truepf(:,5)), std(low_heur_close2truepf(:,5)), mean(low_heur_close2truepf(:,6)), std(low_heur_close2truepf(:,6));
    mean(high_heur_close2truepf(:,1)), std(high_heur_close2truepf(:,1)), mean(high_heur_close2truepf(:,2)), std(high_heur_close2truepf(:,2)), mean(high_heur_close2truepf(:,3)), std(high_heur_close2truepf(:,3)), mean(high_heur_close2truepf(:,4)), std(high_heur_close2truepf(:,4)), mean(high_heur_close2truepf(:,5)), std(high_heur_close2truepf(:,5)), mean(high_heur_close2truepf(:,6)), std(high_heur_close2truepf(:,6))];

%% Functions
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

function heur_intrness_array = compute_heur_intrness_array(heur_thresh, mindist_truepf_thresh, c22thresh, volfracthresh, feasthresh, connthresh, stiffratthresh)
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

function [objs_pen_combined, objs_true_combined, constraints_combined, heuristics_combined, designs_combined] = obtain_combined_data_case(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_feas_bools, constrad_prob_read, repeat3x3_bool, sidenodenum, n_pop, n_runs)
    
    n_members_total = nchoosek(sidenodenum^2,2);    
    data_array = zeros(n_pop,10,n_runs);
    if constrad_prob_read
        designs_array = strings(n_pop,n_runs);
        for i = 1:n_runs
            [data_array(:,:,i), designs_array(:,i)] = read_csv_data(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_feas_bools, constrad_prob_read, repeat3x3_bool, n_members_total, i-1);
        end
    else
        designs_array = zeros(n_pop,n_members_total,n_runs);
        for i = 1:n_runs
            [data_array(:,:,i), designs_array(:,:,i)] = read_csv_data(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_feas_bools, constrad_prob_read, repeat3x3_bool, n_members_total, i-1);
        end
    end
    [pen_obj1_combined, pen_obj2_combined, true_obj1_combined, true_obj2_combined, constraints_combined, heuristics_combined, designs_combined] = create_combined_arrays(data_array, designs_array, constrad_prob_read, n_members_total, n_runs);
    objs_pen_combined = [pen_obj1_combined, pen_obj2_combined];
    objs_true_combined = [true_obj1_combined, true_obj2_combined];
    
end

function [data_array, design_array] = read_csv_data(fib_stiff, partcoll_bools, nodalprop_bools, orient_bools, feas_bools, constrad_read, repeat3x3_bool, n_total_members, run_num)
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

