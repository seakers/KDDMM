%% Observe changes in heuristic indices as a function of termination NFE for GAs
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
sidenum = 3;
nucFac = 3;
biasFactor = 1;
collapsibilityBiasFac = 0.5;
choice_of_model = "Truss"; % "Fibre" -> fibre model, "Truss" -> truss model, "Beam" -> beam model

n_total_members = nchoosek(sidenum^2,2);
n_repeated_members = 2*nchoosek(sidenum,2);
n_variables = n_total_members - n_repeated_members;

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

%% Generate random architectures, objectives, constraints and heuristics to append to GA data
n_des_rand = 100; % number of architectures to generate
n_runs = 10; % number of runs to generate "n_des" architectures

% random_mode = "PartcollAndConn";
% if ~truss_problem
%     random_mode = "Conn";
% end
random_mode = "Random";

add_random_data = true; % true -> add 100 randomly generated designs to the dataset each trial, false -> only use GA data

des_rand_allruns = struct;
objs_rand_allruns = struct;
objs_norm_rand_allruns = struct;
constr_rand_allruns = struct;
heur_rand_allruns = struct;

for i = 1:n_runs
    if add_random_data
        [bool_des_rand, obj_rand, obj_norm_rand, constr_rand, heur_rand] = generate_biased_random_population(n_des_rand, choice_of_model, truss_problem, random_mode, sidenum, sel, r, E, c_ratio, collapsibilityBiasFac, biasFactor);
    else
        bool_des_rand = [];
        obj_rand = [];
        obj_norm_rand = [];
        constr_rand = [];
        heur_rand = [];
    end

    current_field = strcat('trial',num2str(i));
    des_rand_allruns.(current_field) = bool_des_rand;
    objs_rand_allruns.(current_field) = obj_rand;
    objs_norm_rand_allruns.(current_field) = obj_norm_rand;
    constr_rand_allruns.(current_field) = constr_rand;
    heur_rand_allruns.(current_field) = heur_rand;
end

%% Compute heuristic indices for different termination indices for the GA runs
termination_nfes = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000];
I_partcoll = zeros(size(termination_nfes, 2), n_runs);
I_nodalprop = zeros(size(termination_nfes, 2), n_runs);
I_orient = zeros(size(termination_nfes, 2), n_runs);
I_inters = zeros(size(termination_nfes, 2), n_runs);

for i = 1:size(termination_nfes, 2)
    I_heurs_i = compute_heuristic_indices(truss_problem, choice_of_model, add_random_data, termination_nfes(i), case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_inters_bools, constrad_read, sidenum, n_runs, objs_norm_rand_allruns, objs_rand_allruns, constr_rand_allruns, heur_rand_allruns);
    I_partcoll(i,:) = I_heurs_i(:,1);
    I_nodalprop(i,:) = I_heurs_i(:,2);
    I_orient(i,:) = I_heurs_i(:,3);
    I_inters(i,:) = I_heurs_i(:,4);
end

% PLotting
figure
subplot(2,2,1)
errorbar(termination_nfes, mean(I_partcoll,2,'omitnan')', std(I_partcoll,0,2,'omitnan')')
hold on
plot(termination_nfes, zeros(1, length(termination_nfes)))
hold off
ylim([-0.7 0.7])
xlabel("Termination NFEs for GA runs")
ylabel("Heuristic Index")
title("Partial Collapsibility")
subplot(2,2,2)
errorbar(termination_nfes, mean(I_nodalprop,2,'omitnan')', std(I_nodalprop,0,2,'omitnan')')
hold on
plot([0, termination_nfes], zeros(1, (length(termination_nfes) +1)))
hold off
ylim([-0.7 0.7])
xlabel("Termination NFEs for GA runs")
ylabel("Heuristic Index")
title("Nodal Properties")
subplot(2,2,3)
errorbar(termination_nfes, mean(I_orient,2,'omitnan')', std(I_orient,0,2,'omitnan')')
hold on
plot([0, termination_nfes], zeros(1, (length(termination_nfes) +1)))
hold off
ylim([-0.7 0.7])
xlabel("Termination NFEs for GA runs")
ylabel("Heuristic Index")
title("Orientation")
subplot(2,2,4)
errorbar(termination_nfes, mean(I_inters,2,'omitnan')', std(I_inters,0,2,'omitnan')')
hold on
plot([0, termination_nfes], zeros(1, (length(termination_nfes) +1)))
hold off
ylim([-0.7 0.7])
xlabel("Termination NFEs for GA runs")
ylabel("Heuristic Index")
title("Intersection")

%% Functions

function [I_heurs] = compute_heuristic_indices(prob_truss, model_choice, use_random_data, termination_nfe, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, read_constrad, sidenodenum, num_runs, f_norm_rand_allruns, f_true_rand_allruns, constrs_rand_allruns, heurs_rand_allruns)
    % I_heurs = [I_partcoll; I_nodalprop; I_orient; I_inters] * num_rums
    
    % Extract and store GA data upto the termination NFE
    if termination_nfe ~= 0
        pop_size = 100;
        [f_norm_nonans_allgas, f_true_nonans_gas, constr_nonans_allgas, heur_nonans_allgas, ~] = obtain_combined_ga_data_allruns(prob_truss, model_choice, termination_nfe, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, read_constrad, sidenodenum, pop_size, num_runs);   
    end
    
    % Compute heuristic indices for the given termination NFE
    I_partcoll_allruns = zeros(num_runs, 1);
    I_nodalprop_allruns = zeros(num_runs, 1);
    I_orient_allruns = zeros(num_runs, 1);
    I_inters_allruns = zeros(num_runs, 1);

    if termination_nfe ~= 0
        for i = 1:num_runs
            current_field = strcat('trial',num2str(i));
            
            f_norm_nonans_currentcase = f_norm_nonans_allgas.(current_field);
            f_true_nonans_currentcase = f_true_nonans_gas.(current_field);
            constr_nonans_currentcase = constr_nonans_allgas.(current_field);
            heur_nonans_currentcase = heur_nonans_allgas.(current_field);
            
            if use_random_data
                f_norm_rand_trial = f_norm_rand_allruns.(current_field);
                f_rand_trial = f_true_rand_allruns.(current_field);
                constr_rand_trial = constrs_rand_allruns.(current_field);
                heur_rand_trial = heurs_rand_allruns.(current_field);
                
                % Normalized objectives
                obj1_norm_total = cat(1,f_norm_rand_trial(:,1),f_norm_nonans_currentcase(:,1));
                obj2_norm_total = cat(1,f_norm_rand_trial(:,2),f_norm_nonans_currentcase(:,2));
                
                % True objectives
                obj1_total = cat(1,f_rand_trial(:,1),f_true_nonans_currentcase(:,1));
                obj2_total = cat(1,f_rand_trial(:,2),f_true_nonans_currentcase(:,2));
                
                % Constraints
                feas_total = cat(1,constr_rand_trial(:,1),constr_nonans_currentcase(:,1));
                conn_total = cat(1,constr_rand_trial(:,2),constr_nonans_currentcase(:,2));
                if prob_truss
                    stiffrat_total = cat(1,constr_rand_trial(:,3),constr_nonans_currentcase(:,3));
                end
                
                % Heuristics
                coll_total = cat(1,heur_rand_trial(:,1),heur_nonans_currentcase(:,1));
                nod_total = cat(1,heur_rand_trial(:,2),heur_nonans_currentcase(:,2));
                orient_total = cat(1,heur_rand_trial(:,3),heur_nonans_currentcase(:,3));
                %inters_total = feas_total;
                inters_total = cat(1,heur_rand_trial(:,4),heur_nonans_currentcase(:,4));
                
            else
                % Normalized objectives
                obj1_norm_total = f_norm_nonans_currentcase(:,1);
                obj2_norm_total = f_norm_nonans_currentcase(:,2);
                
                % True objectives
                obj1_total = f_true_nonans_currentcase(:,1);
                obj2_total = f_true_nonans_currentcase(:,2);
                
                % Constraints
                feas_total = constr_nonans_currentcase(:,1);
                conn_total = constr_nonans_currentcase(:,2);
                if prob_truss
                    stiffrat_total = constr_nonans_currentcase(:,3);
                end
                
                % Heuristics
                coll_total = heur_nonans_currentcase(:,1);
                nod_total = heur_nonans_currentcase(:,2);
                orient_total = heur_nonans_currentcase(:,3);
                inters_total = heur_nonans_currentcase(:,4);
            end
            
            % Compute normalized and true Pareto Fronts
            norm_objs_pareto = compute_pareto_front(obj1_norm_total,obj2_norm_total);
            true_objs_pareto = compute_pareto_front(-obj1_total,obj2_total);
            true_objs_pareto_correct = [-true_objs_pareto(:,1),true_objs_pareto(:,2)];
            
            fracsat_normpf = size(norm_objs_pareto, 1)/size(obj1_total,1);
            fracsat_feas = length(feas_total(feas_total==1))/size(obj1_total,1);
            fracsat_feas = max([1e-3, fracsat_feas]);
            fracsat_conn = length(feas_total(conn_total==1))/size(obj1_total,1);
            if prob_truss
                fracsat_stiffrat = length(stiffrat_total(stiffrat_total==0))/size(obj1_total,1);
                fracsat_stiffrat = max([1e-3, fracsat_stiffrat]);
            end
            
            % Normalizing objectives and pfs wrt max and min from objectives
            obj1_max = max(obj1_total);
            obj1_min = min(obj1_total);
            obj1_norm_max = max(obj1_norm_total);
            obj1_norm_min = min(obj1_norm_total);
            
            obj2_max = max(obj2_total);
            obj2_min = min(obj2_total);
            obj2_norm_max = max(obj2_norm_total);
            obj2_norm_min = min(obj2_norm_total);
            
            obj1_norm_normalised_total = (obj1_norm_total - obj1_norm_min)/(obj1_norm_max - obj1_norm_min);
            obj2_norm_normalised_total = (obj2_norm_total - obj2_norm_min)/(obj2_norm_max - obj2_norm_min);
            
            obj1_normalised_total = (obj1_total - obj1_min)/(obj1_max - obj1_min);
            obj2_normalised_total = (obj2_total - obj2_min)/(obj2_max - obj2_min);
            
            norm_objs_normalised_pareto = [(norm_objs_pareto(:,1) - obj1_norm_min)/(obj1_norm_max - obj1_norm_min), (norm_objs_pareto(:,2) - obj2_norm_min)/(obj2_norm_max - obj2_norm_min)];
            true_objs_normalised_pareto = [(true_objs_pareto_correct(:,1) - obj1_min)/(obj1_max - obj1_min), (true_objs_pareto_correct(:,2) - obj2_min)/(obj2_max - obj2_min)];
            
            % Compute minimum distance to normalized and true Pareto Fronts
            min_dist_norm_pf_total = zeros(size(obj1_norm_total,1),1);
            min_dist_true_pf_total = zeros(size(obj1_norm_total,1),1);
            for k = 1:size(obj1_norm_total,1)
                min_dist_norm_pf_total(k,1) = compute_min_pf_dist([obj1_norm_normalised_total(k,1),obj2_norm_normalised_total(k,1)],norm_objs_normalised_pareto);
                min_dist_true_pf_total(k,1) = compute_min_pf_dist([obj1_normalised_total(k,1),obj2_normalised_total(k,1)],true_objs_normalised_pareto);
            end
            
            % Computing Pearson's Coefficients
            pearson_coll_normpfdist = corr(coll_total,min_dist_norm_pf_total,'Type','Pearson','Rows','complete');
            pearson_coll_feas = corr(coll_total,feas_total,'Type','Pearson','Rows','complete');
            pearson_coll_conn = corr(coll_total,conn_total,'Type','Pearson','Rows','complete');
            
            pearson_nod_normpfdist = corr(nod_total,min_dist_norm_pf_total,'Type','Pearson','Rows','complete');
            pearson_nod_feas = corr(nod_total,feas_total,'Type','Pearson','Rows','complete');
            pearson_nod_conn = corr(nod_total,conn_total,'Type','Pearson','Rows','complete');
            
            pearson_orient_normpfdist = corr(orient_total,min_dist_norm_pf_total,'Type','Pearson','Rows','complete');
            pearson_orient_feas = corr(orient_total,feas_total,'Type','Pearson','Rows','complete');
            pearson_orient_conn = corr(orient_total,conn_total,'Type','Pearson','Rows','complete');
            
            pearson_inters_normpfdist = corr(inters_total,min_dist_norm_pf_total,'Type','Pearson','Rows','complete');
            pearson_inters_feas = corr(inters_total,feas_total,'Type','Pearson','Rows','complete');
            pearson_inters_conn = corr(inters_total,conn_total,'Type','Pearson','Rows','complete');
            
            % Computing Spearman's Coefficients
            spearman_coll_normpfdist = corr(coll_total,min_dist_norm_pf_total,'Type','Spearman','Rows','complete');
            spearman_coll_feas = corr(coll_total,feas_total,'Type','Spearman','Rows','complete');
            spearman_coll_conn = corr(coll_total,conn_total,'Type','Spearman','Rows','complete');
            
            spearman_nod_normpfdist = corr(nod_total,min_dist_norm_pf_total,'Type','Spearman','Rows','complete');
            spearman_nod_feas = corr(nod_total,feas_total,'Type','Spearman','Rows','complete');
            spearman_nod_conn = corr(nod_total,conn_total,'Type','Spearman','Rows','complete');
            
            spearman_orient_normpfdist = corr(orient_total,min_dist_norm_pf_total,'Type','Spearman','Rows','complete');
            spearman_orient_feas = corr(orient_total,feas_total,'Type','Spearman','Rows','complete');
            spearman_orient_conn = corr(orient_total,conn_total,'Type','Spearman','Rows','complete');
            
            spearman_inters_normpfdist = corr(inters_total,min_dist_norm_pf_total,'Type','Spearman','Rows','complete');
            spearman_inters_feas = corr(inters_total,feas_total,'Type','Spearman','Rows','complete');
            spearman_inters_conn = corr(inters_total,conn_total,'Type','Spearman','Rows','complete');
            
            if prob_truss
                pearson_coll_stiffrat = corr((coll_total),stiffrat_total,'Type','Pearson','Rows','complete');
                pearson_nod_stiffrat = corr(nod_total,stiffrat_total,'Type','Pearson','Rows','complete');
                pearson_orient_stiffrat = corr(orient_total,stiffrat_total,'Type','Pearson','Rows','complete');
                pearson_inters_stiffrat = corr(inters_total,stiffrat_total,'Type','Pearson','Rows','complete');
                spearman_coll_stiffrat = corr((coll_total),stiffrat_total,'Type','Spearman','Rows','complete');
                spearman_nod_stiffrat = corr(nod_total,stiffrat_total,'Type','Spearman','Rows','complete');
                spearman_orient_stiffrat = corr(orient_total,stiffrat_total,'Type','Spearman','Rows','complete');
                spearman_inters_stiffrat = corr(inters_total,stiffrat_total,'Type','Spearman','Rows','complete');
            end
            
            % Calculation of average correlation coefficients for Heuristic index calculations
            corr_avg_coll_normpf = (pearson_coll_normpfdist + spearman_coll_normpfdist)/2;
            corr_avg_coll_feas = (pearson_coll_feas + spearman_coll_feas)/2;
            corr_avg_coll_conn = (pearson_coll_conn + spearman_coll_conn)/2;
            
            corr_avg_nod_normpf = (pearson_nod_normpfdist + spearman_nod_normpfdist)/2;
            corr_avg_nod_feas = (pearson_nod_feas + spearman_nod_feas)/2;
            corr_avg_nod_conn = (pearson_nod_conn + spearman_nod_conn)/2;
            
            corr_avg_orient_normpf = (pearson_orient_normpfdist + spearman_orient_normpfdist)/2;
            corr_avg_orient_feas = (pearson_orient_feas + spearman_orient_feas)/2;
            corr_avg_orient_conn = (pearson_orient_conn + spearman_orient_conn)/2;
            
            corr_avg_inters_normpf = (pearson_inters_normpfdist + spearman_inters_normpfdist)/2;
            corr_avg_inters_feas = (pearson_inters_feas + spearman_inters_feas)/2;
            corr_avg_inters_conn = (pearson_inters_conn + spearman_inters_conn)/2;
            
            if prob_truss
                corr_avg_coll_stiffrat = (pearson_coll_stiffrat + spearman_coll_stiffrat)/2;
                corr_avg_nod_stiffrat = (pearson_nod_stiffrat + spearman_nod_stiffrat)/2;
                corr_avg_orient_stiffrat = (pearson_orient_stiffrat + spearman_orient_stiffrat)/2;
                corr_avg_inters_stiffrat = (pearson_inters_stiffrat + spearman_inters_stiffrat)/2;
            end
            
            if prob_truss
                corr_exp_partcoll = [-1,1,1,-1]; % [normpfdist, feas, conn, stiffrat]
                corr_exp_nodalprop = [-1,1,1,-1]; % [normpfdist, feas, conn, stiffrat]
                corr_exp_orient = [-1,1,1,-1]; % [normpfdist, feas, conn, stiffrat]
                corr_exp_inters = [-1,1,1,-1]; % [normpfdist, feas, conn, stiffrat]
            else
                corr_exp_partcoll = [-1,1,1]; % [normpfdist, feas, conn]
                corr_exp_nodalprop = [-1,1,1]; % [normpfdist, feas, conn]
                corr_exp_orient = [-1,1,1]; % [normpfdist, feas, conn]
                corr_exp_inters = [-1,1,1]; % [normpfdist, feas, conn]
            end
            
            I_partcoll_allruns(i,1) = compute_heuristic_I1_contribution_run(corr_avg_coll_normpf, corr_exp_partcoll(1), fracsat_normpf) + ...
                compute_heuristic_I1_contribution_run(corr_avg_coll_feas, corr_exp_partcoll(2), fracsat_feas) + ...
                compute_heuristic_I1_contribution_run(corr_avg_coll_conn, corr_exp_partcoll(3), fracsat_conn);
            I_nodalprop_allruns(i,1) = compute_heuristic_I1_contribution_run(corr_avg_nod_normpf, corr_exp_nodalprop(1), fracsat_normpf) + ...
                compute_heuristic_I1_contribution_run(corr_avg_nod_feas, corr_exp_nodalprop(2), fracsat_feas) + ...
                compute_heuristic_I1_contribution_run(corr_avg_nod_conn, corr_exp_nodalprop(3), fracsat_conn);
            I_orient_allruns(i,1) = compute_heuristic_I1_contribution_run(corr_avg_orient_normpf, corr_exp_orient(1), fracsat_normpf) + ...
                compute_heuristic_I1_contribution_run(corr_avg_orient_feas, corr_exp_orient(2), fracsat_feas) + ...
                compute_heuristic_I1_contribution_run(corr_avg_orient_conn, corr_exp_orient(3), fracsat_conn);
            I_inters_allruns(i,1) = compute_heuristic_I1_contribution_run(corr_avg_inters_normpf, corr_exp_inters(1), fracsat_normpf) + ...
                compute_heuristic_I1_contribution_run(corr_avg_inters_feas, corr_exp_inters(2), fracsat_feas) + ...
                compute_heuristic_I1_contribution_run(corr_avg_inters_conn, corr_exp_inters(3), fracsat_conn);
            
            if prob_truss
                I_partcoll_allruns(i,1) = I_partcoll_allruns(i,1) + compute_heuristic_I1_contribution_run(corr_avg_coll_stiffrat, corr_exp_partcoll(4), fracsat_stiffrat);
                I_nodalprop_allruns(i,1) = I_nodalprop_allruns(i,1) + compute_heuristic_I1_contribution_run(corr_avg_nod_stiffrat, corr_exp_nodalprop(4), fracsat_stiffrat);
                I_orient_allruns(i,1) = I_orient_allruns(i,1) + compute_heuristic_I1_contribution_run(corr_avg_orient_stiffrat, corr_exp_orient(4), fracsat_stiffrat);
                I_inters_allruns(i,1) = I_inters_allruns(i,1) + compute_heuristic_I1_contribution_run(corr_avg_inters_stiffrat, corr_exp_inters(4), fracsat_stiffrat);
                
                I_partcoll_allruns(i,1) = I_partcoll_allruns(i,1)/4;
                I_nodalprop_allruns(i,1) = I_nodalprop_allruns(i,1)/4;
                I_orient_allruns(i,1) = I_orient_allruns(i,1)/4;
                I_inters_allruns(i,1) = I_inters_allruns(i,1)/4;
            else
                I_partcoll_allruns(i,1) = I_partcoll_allruns(i,1)/3;
                I_nodalprop_allruns(i,1) = I_nodalprop_allruns(i,1)/3;
                I_orient_allruns(i,1) = I_orient_allruns(i,1)/3;
                I_inters_allruns(i,1) = I_inters_allruns(i,1)/3;
            end
        end
    end
    
    I_heurs = [I_partcoll_allruns, I_nodalprop_allruns, I_orient_allruns, I_inters_allruns];
end

function [objs_pen_nonans_allcases, objs_true_nonans_allcases, constraints_nonans_allcases, heuristics_nonans_allcases, designs_nonans_allcases] = obtain_combined_ga_data_allruns(truss_prob, model_used, term_nfe, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, constrad_prob_read, sidenodenum, n_pop, n_runs)
    
    n_members_total = nchoosek(sidenodenum^2,2);    
   
    objs_pen_nonans_allcases = struct;
    objs_true_nonans_allcases = struct;
    constraints_nonans_allcases = struct;
    heuristics_nonans_allcases = struct;
    designs_nonans_allcases = struct;
    
    for i = 1:n_runs
        [data_array, designs_array] = read_csv_data_tillnfe(truss_prob, model_used, constrad_prob_read, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_inters_bools, term_nfe, n_members_total, i-1);
        
        n_constr = size(data_array,2) - 1 - 4 - 4; % number of constraints changes based on problem, so subtract 4 (2 pen. and 2 true objs.) and 4 (heurs) from number of total columns, + 1 is for the NFE column
        data_array_nonans_bool = any(isnan(data_array),2);
        data_array_nonans = data_array(~data_array_nonans_bool,:);
        designs_array_nonans = designs_array(~data_array_nonans_bool,:);
        
        current_field = strcat('trial',num2str(i));
        
        objs_pen_nonans_allcases.(current_field) = [data_array_nonans(:,2), data_array_nonans(:,3)];
        objs_true_nonans_allcases.(current_field) = [data_array_nonans(:,4), data_array_nonans(:,5)];
        constraints_nonans_allcases.(current_field) = data_array_nonans(:,6:6+n_constr-1);
        heuristics_nonans_allcases.(current_field) = data_array_nonans(:,6+n_constr:end);
        designs_nonans_allcases.(current_field) = designs_array_nonans;
        
    end
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
        filepath_moea = "Epsilon MOEA - Metrics\\";
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

function index_contribution = compute_heuristic_I1_contribution2(corr_array_heur_param, idx_corr_heur_param, supp_array_param)
    % corr_array_heur_param and supp_array_param are (n x 1) arrays where n is the number of runs
    log_arg = idx_corr_heur_param*mean(corr_array_heur_param);
    index_contribution = log_arg*(-1*log10(mean(supp_array_param)));
end

function index_contribution = compute_heuristic_I1_contribution_run(corr_heur_param, idx_corr_heur_param, supp_param)
    % corr_heur_param and supp_param are scalars
    log_arg = idx_corr_heur_param*corr_heur_param;
    index_contribution = log_arg*(-1*log10(supp_param));
end