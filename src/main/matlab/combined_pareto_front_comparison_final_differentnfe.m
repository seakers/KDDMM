%% Plot combined pareto front comparison
clear all
close all
clc

%% CSV Read parameters
fibre_model_used = false;
sidenum = 3;

num_runs = 30; % change based on run 
pop_size = 100; % change based on run

E = 1.8162e6; % Young's Modulus for polymeric material (example: 1.8162 MPa for SIL material)

% partcoll_bools = [int_pen, AOS, biased_init, ACH] boolean array
% nodalprop_bools = [int_pen, AOS, biased_init, ACH] boolean array
% orient_bools = [int_pen, AOS, biased_init, ACH] boolean array
% feas_bools = [biased_init, AOS] boolean array

%% Comparing Combined Pareto Fronts (eps-MOEA vs AOS - Orientation)(constant radii problem)
% Case 1 - Epsilon MOEA
constrad_read_case1 = true;
repeat3x3_case1 = true;
case1_partcoll_bools = [false, false, false, false];
case1_nodalprop_bools = [false, false, false, false];
case1_orient_bools = [false, false, false, false];
case1_feas_bools = [false, false];

% Case 2 - AOS - Orientation + Intersection, Bias Init - Intersection
constrad_read_case2 = true;
repeat3x3_case2 = true;
case2_partcoll_bools = [false, false, false, false];
case2_nodalprop_bools = [false, false, false, false];
case2_orient_bools = [false, true, false, false];
case2_feas_bools = [true, true];

% Case 3 - AOS - Orientation 
constrad_read_case3 = true;
repeat3x3_case3 = true;
case3_partcoll_bools = [false, false, false, false];
case3_nodalprop_bools = [false, false, false, false];
case3_orient_bools = [false, true, false, false];
case3_feas_bools = [false, false];

% Case 4 - AOS - Partial Collapsibility and Nodal Properties
constrad_read_case4 = true;
repeat3x3_case4 = true;
case4_partcoll_bools = [false, true, false, false];
case4_nodalprop_bools = [false, true, false, false];
case4_orient_bools = [false, false, false, false];
case4_feas_bools = [false, false];

% Generate combined true and penalized pareto front arrays for each run
% case for different nfe thresholds

% NFE threshold = 250
nfe_thresh1 = 250;
% Case 1
[f_pen_pareto_combined_case1_thresh1, f_true_pareto_combined_case1_thresh1, f_true_feas_pareto_combined_case1_thresh1, des_pareto_combined_case1_thresh1] = obtain_combined_pareto_data_case(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_feas_bools, constrad_read_case1, repeat3x3_case1, nfe_thresh1, sidenum, pop_size, num_runs); 

% Case 2
[f_pen_pareto_combined_case2_thresh1, f_true_pareto_combined_case2_thresh1, f_true_feas_pareto_combined_case2_thresh1, des_pareto_combined_case2_thresh1] = obtain_combined_pareto_data_case(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_feas_bools, constrad_read_case2, repeat3x3_case2, nfe_thresh1, sidenum, pop_size, num_runs); 

% Case 3
[f_pen_pareto_combined_case3_thresh1, f_true_pareto_combined_case3_thresh1, f_true_feas_pareto_combined_case3_thresh1, des_pareto_combined_case3_thresh1] = obtain_combined_pareto_data_case(fibre_model_used, case3_partcoll_bools, case3_nodalprop_bools, case3_orient_bools, case3_feas_bools, constrad_read_case3, repeat3x3_case3, nfe_thresh1, sidenum, pop_size, num_runs); 

% Case 4
[f_pen_pareto_combined_case4_thresh1, f_true_pareto_combined_case4_thresh1, f_true_feas_pareto_combined_case4_thresh1, des_pareto_combined_case4_thresh1] = obtain_combined_pareto_data_case(fibre_model_used, case4_partcoll_bools, case4_nodalprop_bools, case4_orient_bools, case4_feas_bools, constrad_read_case4, repeat3x3_case4, nfe_thresh1, sidenum, pop_size, num_runs); 

penobj1_min_case1_thresh1 = min(f_pen_pareto_combined_case1_thresh1(:,1));
penobj1_min_case2_thresh1 = min(f_pen_pareto_combined_case2_thresh1(:,1));
penobj1_min_case3_thresh1 = min(f_pen_pareto_combined_case3_thresh1(:,1));
penobj1_min_case4_thresh1 = min(f_pen_pareto_combined_case4_thresh1(:,1));
penobj1_min_thresh1 = min([penobj1_min_case1_thresh1,penobj1_min_case2_thresh1,penobj1_min_case3_thresh1,penobj1_min_case4_thresh1]);

penobj2_min_case1_thresh1 = min(f_pen_pareto_combined_case1_thresh1(:,2));
penobj2_min_case2_thresh1 = min(f_pen_pareto_combined_case2_thresh1(:,2));
penobj2_min_case3_thresh1 = min(f_pen_pareto_combined_case3_thresh1(:,2));
penobj2_min_case4_thresh1 = min(f_pen_pareto_combined_case4_thresh1(:,2));
penobj2_min_thresh1 = min([penobj2_min_case1_thresh1,penobj2_min_case2_thresh1,penobj2_min_case3_thresh1,penobj2_min_case4_thresh1]);

trueobj1_max_case1_thresh1 = max(f_true_pareto_combined_case1_thresh1(:,1));
trueobj1_max_case2_thresh1 = max(f_true_pareto_combined_case2_thresh1(:,1));
trueobj1_max_case3_thresh1 = max(f_true_pareto_combined_case3_thresh1(:,1));
trueobj1_max_case4_thresh1 = max(f_true_pareto_combined_case4_thresh1(:,1));
trueobj1_max_thresh1 = max([trueobj1_max_case1_thresh1,trueobj1_max_case2_thresh1,trueobj1_max_case3_thresh1,trueobj1_max_case4_thresh1]);

trueobj2_min_case1_thresh1 = min(f_true_pareto_combined_case1_thresh1(:,2));
trueobj2_min_case2_thresh1 = min(f_true_pareto_combined_case2_thresh1(:,2));
trueobj2_min_case3_thresh1 = min(f_true_pareto_combined_case3_thresh1(:,2));
trueobj2_min_case4_thresh1 = min(f_true_pareto_combined_case4_thresh1(:,2));
trueobj2_min_thresh1 = min([trueobj2_min_case1_thresh1,trueobj2_min_case2_thresh1,trueobj2_min_case3_thresh1,trueobj2_min_case4_thresh1]);

utopia_pen_thresh1 = [penobj1_min_thresh1, penobj2_min_thresh1];
utopia_true_thresh1 = [trueobj1_max_thresh1, trueobj2_min_thresh1];
    
% NFE threshold = 500
nfe_thresh2 = 500;
% Case 1
[f_pen_pareto_combined_case1_thresh2, f_true_pareto_combined_case1_thresh2, f_true_feas_pareto_combined_case1_thresh2, des_pareto_combined_case1_thresh2] = obtain_combined_pareto_data_case(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_feas_bools, constrad_read_case1, repeat3x3_case1, nfe_thresh2, sidenum, pop_size, num_runs); 

% Case 2
[f_pen_pareto_combined_case2_thresh2, f_true_pareto_combined_case2_thresh2, f_true_feas_pareto_combined_case2_thresh2, des_pareto_combined_case2_thresh2] = obtain_combined_pareto_data_case(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_feas_bools, constrad_read_case2, repeat3x3_case2, nfe_thresh2, sidenum, pop_size, num_runs); 

% Case 3
[f_pen_pareto_combined_case3_thresh2, f_true_pareto_combined_case3_thresh2, f_true_feas_pareto_combined_case3_thresh2, des_pareto_combined_case3_thresh2] = obtain_combined_pareto_data_case(fibre_model_used, case3_partcoll_bools, case3_nodalprop_bools, case3_orient_bools, case3_feas_bools, constrad_read_case3, repeat3x3_case3, nfe_thresh2, sidenum, pop_size, num_runs); 

% Case 4
[f_pen_pareto_combined_case4_thresh2, f_true_pareto_combined_case4_thresh2, f_true_feas_pareto_combined_case4_thresh2, des_pareto_combined_case4_thresh2] = obtain_combined_pareto_data_case(fibre_model_used, case4_partcoll_bools, case4_nodalprop_bools, case4_orient_bools, case4_feas_bools, constrad_read_case4, repeat3x3_case4, nfe_thresh2, sidenum, pop_size, num_runs); 

penobj1_min_case1_thresh2 = min(f_pen_pareto_combined_case1_thresh2(:,1));
penobj1_min_case2_thresh2 = min(f_pen_pareto_combined_case2_thresh2(:,1));
penobj1_min_case3_thresh2 = min(f_pen_pareto_combined_case3_thresh2(:,1));
penobj1_min_case4_thresh2 = min(f_pen_pareto_combined_case4_thresh2(:,1));
penobj1_min_thresh2 = min([penobj1_min_case1_thresh2,penobj1_min_case2_thresh2,penobj1_min_case3_thresh2,penobj1_min_case4_thresh2]);

penobj2_min_case1_thresh2 = min(f_pen_pareto_combined_case1_thresh2(:,2));
penobj2_min_case2_thresh2 = min(f_pen_pareto_combined_case2_thresh2(:,2));
penobj2_min_case3_thresh2 = min(f_pen_pareto_combined_case3_thresh2(:,2));
penobj2_min_case4_thresh2 = min(f_pen_pareto_combined_case4_thresh2(:,2));
penobj2_min_thresh2 = min([penobj2_min_case1_thresh2,penobj2_min_case2_thresh2,penobj2_min_case3_thresh2,penobj2_min_case4_thresh2]);

trueobj1_max_case1_thresh2 = max(f_true_pareto_combined_case1_thresh2(:,1));
trueobj1_max_case2_thresh2 = max(f_true_pareto_combined_case2_thresh2(:,1));
trueobj1_max_case3_thresh2 = max(f_true_pareto_combined_case3_thresh2(:,1));
trueobj1_max_case4_thresh2 = max(f_true_pareto_combined_case4_thresh2(:,1));
trueobj1_max_thresh2 = max([trueobj1_max_case1_thresh2,trueobj1_max_case2_thresh2,trueobj1_max_case3_thresh2,trueobj1_max_case4_thresh2]);

trueobj2_min_case1_thresh2 = min(f_true_pareto_combined_case1_thresh2(:,2));
trueobj2_min_case2_thresh2 = min(f_true_pareto_combined_case2_thresh2(:,2));
trueobj2_min_case3_thresh2 = min(f_true_pareto_combined_case3_thresh2(:,2));
trueobj2_min_case4_thresh2 = min(f_true_pareto_combined_case4_thresh2(:,2));
trueobj2_min_thresh2 = min([trueobj2_min_case1_thresh2,trueobj2_min_case2_thresh2,trueobj2_min_case3_thresh2,trueobj2_min_case4_thresh2]);

utopia_pen_thresh2 = [penobj1_min_thresh2, penobj2_min_thresh2];
utopia_true_thresh2 = [trueobj1_max_thresh2, trueobj2_min_thresh2];

% NFE threshold = 1500
nfe_thresh3 = 1500;
% Case 1
[f_pen_pareto_combined_case1_thresh3, f_true_pareto_combined_case1_thresh3, f_true_feas_pareto_combined_case1_thresh3, des_pareto_combined_case1_thresh3] = obtain_combined_pareto_data_case(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_feas_bools, constrad_read_case1, repeat3x3_case1, nfe_thresh3, sidenum, pop_size, num_runs); 

% Case 2
[f_pen_pareto_combined_case2_thresh3, f_true_pareto_combined_case2_thresh3, f_true_feas_pareto_combined_case2_thresh3, des_pareto_combined_case2_thresh3] = obtain_combined_pareto_data_case(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_feas_bools, constrad_read_case2, repeat3x3_case2, nfe_thresh3, sidenum, pop_size, num_runs); 

% Case 3
[f_pen_pareto_combined_case3_thresh3, f_true_pareto_combined_case3_thresh3, f_true_feas_pareto_combined_case3_thresh3, des_pareto_combined_case3_thresh3] = obtain_combined_pareto_data_case(fibre_model_used, case3_partcoll_bools, case3_nodalprop_bools, case3_orient_bools, case3_feas_bools, constrad_read_case3, repeat3x3_case3, nfe_thresh3, sidenum, pop_size, num_runs); 

% Case 4
[f_pen_pareto_combined_case4_thresh3, f_true_pareto_combined_case4_thresh3, f_true_feas_pareto_combined_case4_thresh3, des_pareto_combined_case4_thresh3] = obtain_combined_pareto_data_case(fibre_model_used, case4_partcoll_bools, case4_nodalprop_bools, case4_orient_bools, case4_feas_bools, constrad_read_case4, repeat3x3_case4, nfe_thresh3, sidenum, pop_size, num_runs); 

penobj1_min_case1_thresh3 = min(f_pen_pareto_combined_case1_thresh3(:,1));
penobj1_min_case2_thresh3 = min(f_pen_pareto_combined_case2_thresh3(:,1));
penobj1_min_case3_thresh3 = min(f_pen_pareto_combined_case3_thresh3(:,1));
penobj1_min_case4_thresh3 = min(f_pen_pareto_combined_case4_thresh3(:,1));
penobj1_min_thresh3 = min([penobj1_min_case1_thresh3,penobj1_min_case2_thresh3,penobj1_min_case3_thresh3,penobj1_min_case4_thresh3]);

penobj2_min_case1_thresh3 = min(f_pen_pareto_combined_case1_thresh3(:,2));
penobj2_min_case2_thresh3 = min(f_pen_pareto_combined_case2_thresh3(:,2));
penobj2_min_case3_thresh3 = min(f_pen_pareto_combined_case3_thresh3(:,2));
penobj2_min_case4_thresh3 = min(f_pen_pareto_combined_case4_thresh3(:,2));
penobj2_min_thresh3 = min([penobj2_min_case1_thresh3,penobj2_min_case2_thresh3,penobj2_min_case3_thresh3,penobj2_min_case4_thresh3]);

trueobj1_max_case1_thresh3 = max(f_true_pareto_combined_case1_thresh3(:,1));
trueobj1_max_case2_thresh3 = max(f_true_pareto_combined_case2_thresh3(:,1));
trueobj1_max_case3_thresh3 = max(f_true_pareto_combined_case3_thresh3(:,1));
trueobj1_max_case4_thresh3 = max(f_true_pareto_combined_case4_thresh3(:,1));
trueobj1_max_thresh3 = max([trueobj1_max_case1_thresh3,trueobj1_max_case2_thresh3,trueobj1_max_case3_thresh3,trueobj1_max_case4_thresh3]);

trueobj2_min_case1_thresh3 = min(f_true_pareto_combined_case1_thresh3(:,2));
trueobj2_min_case2_thresh3 = min(f_true_pareto_combined_case2_thresh3(:,2));
trueobj2_min_case3_thresh3 = min(f_true_pareto_combined_case3_thresh3(:,2));
trueobj2_min_case4_thresh3 = min(f_true_pareto_combined_case4_thresh3(:,2));
trueobj2_min_thresh3 = min([trueobj2_min_case1_thresh3,trueobj2_min_case2_thresh3,trueobj2_min_case3_thresh3,trueobj2_min_case4_thresh3]);

utopia_pen_thresh3 = [penobj1_min_thresh3, penobj2_min_thresh3];
utopia_true_thresh3 = [trueobj1_max_thresh3, trueobj2_min_thresh3];

% NFE threshold = 3000
nfe_thresh4 = 3000;
% Case 1
[f_pen_pareto_combined_case1_thresh4, f_true_pareto_combined_case1_thresh4, f_true_feas_pareto_combined_case1_thresh4, des_pareto_combined_case1_thresh4] = obtain_combined_pareto_data_case(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_feas_bools, constrad_read_case1, repeat3x3_case1, nfe_thresh4, sidenum, pop_size, num_runs); 

% Case 2
[f_pen_pareto_combined_case2_thresh4, f_true_pareto_combined_case2_thresh4, f_true_feas_pareto_combined_case2_thresh4, des_pareto_combined_case2_thresh4] = obtain_combined_pareto_data_case(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_feas_bools, constrad_read_case2, repeat3x3_case2, nfe_thresh4, sidenum, pop_size, num_runs); 

% Case 3
[f_pen_pareto_combined_case3_thresh4, f_true_pareto_combined_case3_thresh4, f_true_feas_pareto_combined_case3_thresh4, des_pareto_combined_case3_thresh4] = obtain_combined_pareto_data_case(fibre_model_used, case3_partcoll_bools, case3_nodalprop_bools, case3_orient_bools, case3_feas_bools, constrad_read_case3, repeat3x3_case3, nfe_thresh4, sidenum, pop_size, num_runs); 

% Case 4
[f_pen_pareto_combined_case4_thresh4, f_true_pareto_combined_case4_thresh4, f_true_feas_pareto_combined_case4_thresh4, des_pareto_combined_case4_thresh4] = obtain_combined_pareto_data_case(fibre_model_used, case4_partcoll_bools, case4_nodalprop_bools, case4_orient_bools, case4_feas_bools, constrad_read_case4, repeat3x3_case4, nfe_thresh4, sidenum, pop_size, num_runs); 

penobj1_min_case1_thresh4 = min(f_pen_pareto_combined_case1_thresh4(:,1));
penobj1_min_case2_thresh4 = min(f_pen_pareto_combined_case2_thresh4(:,1));
penobj1_min_case3_thresh4 = min(f_pen_pareto_combined_case3_thresh4(:,1));
penobj1_min_case4_thresh4 = min(f_pen_pareto_combined_case4_thresh4(:,1));
penobj1_min_thresh4 = min([penobj1_min_case1_thresh4,penobj1_min_case2_thresh4,penobj1_min_case3_thresh4,penobj1_min_case4_thresh4]);

penobj2_min_case1_thresh4 = min(f_pen_pareto_combined_case1_thresh4(:,2));
penobj2_min_case2_thresh4 = min(f_pen_pareto_combined_case2_thresh4(:,2));
penobj2_min_case3_thresh4 = min(f_pen_pareto_combined_case3_thresh4(:,2));
penobj2_min_case4_thresh4 = min(f_pen_pareto_combined_case4_thresh4(:,2));
penobj2_min_thresh4 = min([penobj2_min_case1_thresh4,penobj2_min_case2_thresh4,penobj2_min_case3_thresh4,penobj2_min_case4_thresh4]);

trueobj1_max_case1_thresh4 = max(f_true_pareto_combined_case1_thresh4(:,1));
trueobj1_max_case2_thresh4 = max(f_true_pareto_combined_case2_thresh4(:,1));
trueobj1_max_case3_thresh4 = max(f_true_pareto_combined_case3_thresh4(:,1));
trueobj1_max_case4_thresh4 = max(f_true_pareto_combined_case4_thresh4(:,1));
trueobj1_max_thresh4 = max([trueobj1_max_case1_thresh4,trueobj1_max_case2_thresh4,trueobj1_max_case3_thresh4,trueobj1_max_case4_thresh4]);

trueobj2_min_case1_thresh4 = min(f_true_pareto_combined_case1_thresh4(:,2));
trueobj2_min_case2_thresh4 = min(f_true_pareto_combined_case2_thresh4(:,2));
trueobj2_min_case3_thresh4 = min(f_true_pareto_combined_case3_thresh4(:,2));
trueobj2_min_case4_thresh4 = min(f_true_pareto_combined_case4_thresh4(:,2));
trueobj2_min_thresh4 = min([trueobj2_min_case1_thresh4,trueobj2_min_case2_thresh4,trueobj2_min_case3_thresh4,trueobj2_min_case4_thresh4]);

utopia_pen_thresh4 = [penobj1_min_thresh4, penobj2_min_thresh4];
utopia_true_thresh4 = [trueobj1_max_thresh4, trueobj2_min_thresh4];

% Plotting
figure
subplot(2,2,1)
scatter(f_pen_pareto_combined_case1_thresh1(:,1), f_pen_pareto_combined_case1_thresh1(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_pen_pareto_combined_case2_thresh1(:,1), f_pen_pareto_combined_case2_thresh1(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_pen_pareto_combined_case3_thresh1(:,1), f_pen_pareto_combined_case3_thresh1(:,2), 'Marker', 's', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_pen_pareto_combined_case4_thresh1(:,1), f_pen_pareto_combined_case4_thresh1(:,2), 'Marker', '^', 'MarkerEdgeColor', 'black') 
% hold on
% scatter(utopia_pen_thresh1(1), utopia_pen_thresh1(2), 'Marker', 'p', 'MarkerEdgeColor', 'red')
hold off
xlabel('Penalized $C_{22}$','Interpreter','Latex','FontSize',14)
ylabel('Penalized $v_f$','Interpreter','Latex','FontSize',14)
% xlim([-1,1])
% ylim([0,1])
%legend('No Heur.','Orient+Inters','Orient','PartColl+NodalProp','Location','northeast')
title('250 NFE')

subplot(2,2,2)
scatter(f_pen_pareto_combined_case1_thresh2(:,1), f_pen_pareto_combined_case1_thresh2(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_pen_pareto_combined_case2_thresh2(:,1), f_pen_pareto_combined_case2_thresh2(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_pen_pareto_combined_case3_thresh2(:,1), f_pen_pareto_combined_case3_thresh2(:,2), 'Marker', 's', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_pen_pareto_combined_case4_thresh2(:,1), f_pen_pareto_combined_case4_thresh2(:,2), 'Marker', '^', 'MarkerEdgeColor', 'black') 
% hold on
% scatter(utopia_pen_thresh2(1), utopia_pen_thresh2(2), 'Marker', 'p', 'MarkerEdgeColor', 'red')
hold off
xlabel('Penalized $C_{22}$','Interpreter','Latex','FontSize',14)
ylabel('Penalized $v_f$','Interpreter','Latex','FontSize',14)
% xlim([-1,1])
% ylim([0,1])
%legend('Without Heuristics','With Heuristics','Location','best')
title('500 NFE')

subplot(2,2,3)
scatter(f_pen_pareto_combined_case1_thresh3(:,1), f_pen_pareto_combined_case1_thresh3(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_pen_pareto_combined_case2_thresh3(:,1), f_pen_pareto_combined_case2_thresh3(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_pen_pareto_combined_case3_thresh3(:,1), f_pen_pareto_combined_case3_thresh3(:,2), 'Marker', 's', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_pen_pareto_combined_case4_thresh3(:,1), f_pen_pareto_combined_case4_thresh3(:,2), 'Marker', '^', 'MarkerEdgeColor', 'black') 
% hold on
% scatter(utopia_pen_thresh3(1), utopia_pen_thresh3(2), 'Marker', 'p', 'MarkerEdgeColor', 'red')
hold off
xlabel('Penalized $C_{22}$','Interpreter','Latex','FontSize',14)
ylabel('Penalized $v_f$','Interpreter','Latex','FontSize',14)
% xlim([-1,1])
% ylim([0,1])
%legend('Without Heuristics','With Heuristics','Location','best')
title('1500 NFE')

subplot(2,2,4)
scatter(f_pen_pareto_combined_case1_thresh4(:,1), f_pen_pareto_combined_case1_thresh4(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_pen_pareto_combined_case2_thresh4(:,1), f_pen_pareto_combined_case2_thresh4(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_pen_pareto_combined_case3_thresh4(:,1), f_pen_pareto_combined_case3_thresh4(:,2), 'Marker', 's', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_pen_pareto_combined_case4_thresh4(:,1), f_pen_pareto_combined_case4_thresh4(:,2), 'Marker', '^', 'MarkerEdgeColor', 'black') 
% hold on
% scatter(utopia_pen_thresh4(1), utopia_pen_thresh4(2), 'Marker', 'p', 'MarkerEdgeColor', 'red')
hold off
xlabel('Penalized $C_{22}$','Interpreter','Latex','FontSize',14)
ylabel('Penalized $v_f$','Interpreter','Latex','FontSize',14)
% xlim([-1,1])
% ylim([0,1])
%legend('Without Heuristics','With Heuristics','Location','best')
title('3000 NFE')
saveas(gcf,'pareto_pen_nfequadchart.png')

% Plotting feasible designs in true objectives space
figure
subplot(2,2,1)
scatter(f_true_feas_pareto_combined_case1_thresh1(:,1), f_true_feas_pareto_combined_case1_thresh1(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case2_thresh1(:,1), f_true_feas_pareto_combined_case2_thresh1(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case3_thresh1(:,1), f_true_feas_pareto_combined_case3_thresh1(:,2), 'Marker', 's', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case4_thresh1(:,1), f_true_feas_pareto_combined_case4_thresh1(:,2), 'Marker', '^', 'MarkerEdgeColor', 'black') 
% hold on
% scatter(utopia_true_thresh1(1), utopia_true_thresh1(2), 'Marker', 'p', 'MarkerEdgeColor', 'red')
hold off
ax = gca;
ax.FontSize = 13;
xlabel('$C_{22}$','Interpreter','Latex','FontSize',13)
ylabel('$v_f$','Interpreter','Latex','FontSize',13)
% xlim([0,E/2])
% ylim([0,1])1
%legend('No Heur.','Orient+Inters','Orient','PartColl+NodalProp','Location','northeast')
title('250 NFE','FontSize',13)

subplot(2,2,2)
scatter(f_true_feas_pareto_combined_case1_thresh2(:,1), f_true_feas_pareto_combined_case1_thresh2(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case2_thresh2(:,1), f_true_feas_pareto_combined_case2_thresh2(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case3_thresh2(:,1), f_true_feas_pareto_combined_case3_thresh2(:,2), 'Marker', 's', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case4_thresh2(:,1), f_true_feas_pareto_combined_case4_thresh2(:,2), 'Marker', '^', 'MarkerEdgeColor', 'black') 
% hold on
% scatter(utopia_true_thresh2(1), utopia_true_thresh2(2), 'Marker', 'p', 'MarkerEdgeColor', 'red')
hold off
ax = gca;
ax.FontSize = 13;
xlabel('$C_{22}$','Interpreter','Latex','FontSize',13)
ylabel('$v_f$','Interpreter','Latex','FontSize',13)
% xlim([0,E/2])
% ylim([0,1])
%legend('Without Heuristics','With Heuristics','Location','best')
title('500 NFE','FontSize',13)

subplot(2,2,3)
scatter(f_true_feas_pareto_combined_case1_thresh3(:,1), f_true_feas_pareto_combined_case1_thresh3(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case2_thresh3(:,1), f_true_feas_pareto_combined_case2_thresh3(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case3_thresh3(:,1), f_true_feas_pareto_combined_case3_thresh3(:,2), 'Marker', 's', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case4_thresh3(:,1), f_true_feas_pareto_combined_case4_thresh3(:,2), 'Marker', '^', 'MarkerEdgeColor', 'black') 
% hold on
% scatter(utopia_true_thresh3(1), utopia_true_thresh3(2), 'Marker', 'p', 'MarkerEdgeColor', 'red')
hold off
ax = gca;
ax.FontSize = 13;
xlabel('$C_{22}$','Interpreter','Latex','FontSize',13)
ylabel('$v_f$','Interpreter','Latex','FontSize',13)
% xlim([0,E/2])
% ylim([0,1])
%legend('Without Heuristics','With Heuristics','Location','best')
title('1500 NFE','FontSize',13)

subplot(2,2,4)
scatter(f_true_feas_pareto_combined_case1_thresh4(:,1), f_true_feas_pareto_combined_case1_thresh4(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case2_thresh4(:,1), f_true_feas_pareto_combined_case2_thresh4(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case3_thresh4(:,1), f_true_feas_pareto_combined_case3_thresh4(:,2), 'Marker', 's', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case4_thresh4(:,1), f_true_feas_pareto_combined_case4_thresh4(:,2), 'Marker', '^', 'MarkerEdgeColor', 'black') 
% hold on
% scatter(utopia_true_thresh4(1), utopia_true_thresh4(2), 'Marker', 'p', 'MarkerEdgeColor', 'red')
hold off
ax = gca;
ax.FontSize = 13;
xlabel('$C_{22}$','Interpreter','Latex','FontSize',13)
ylabel('$v_f$','Interpreter','Latex','FontSize',13)
% xlim([0,E/2])
% ylim([0,1])
%legend('Without Heuristics','With Heuristics','Location','best')
title('3000 NFE','FontSize',13)
saveas(gcf,'pareto_truefeas_nfequadchart.png')

%% Comparing Combined Pareto Fronts (eps-MOEA vs AOS - Orientation)(variable radii problem)
% Case 1 - Epsilon MOEA
constrad_read_case1 = false;
repeat3x3_case1 = true;
case1_partcoll_bools = [false, false, false, false];
case1_nodalprop_bools = [false, false, false, false];
case1_orient_bools = [false, false, false, false];
case1_feas_bools = [false, false];

% Case 2 - AOS - PartColl & Orientation, Bias Init - Orientation, ACH -
% NodalProp
constrad_read_case2 = false;
repeat3x3_case2 = true;
case2_partcoll_bools = [false, false, false, false];
case2_nodalprop_bools = [false, false, false, false];
case2_orient_bools = [false, true, false, false];
case2_feas_bools = [true, true];

% Generate combined true and penalized pareto front arrays for each run
% case for different nfe thresholds

% NFE threshold = 250
nfe_thresh1 = 250;
% Case 1
[f_pen_pareto_combined_case1_thresh1, f_true_pareto_combined_case1_thresh1, f_true_feas_pareto_combined_case1_thresh1, des_pareto_combined_case1_thresh1] = obtain_combined_pareto_data_case(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_feas_bools, constrad_read_case1, repeat3x3_case1, nfe_thresh1, sidenum, pop_size, num_runs); 

% Case 2
[f_pen_pareto_combined_case2_thresh1, f_true_pareto_combined_case2_thresh1, f_true_feas_pareto_combined_case2_thresh1, des_pareto_combined_case2_thresh1] = obtain_combined_pareto_data_case(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_feas_bools, constrad_read_case2, repeat3x3_case2, nfe_thresh1, sidenum, pop_size, num_runs); 

% NFE threshold = 500
nfe_thresh2 = 500;
% Case 1
[f_pen_pareto_combined_case1_thresh2, f_true_pareto_combined_case1_thresh2, f_true_feas_pareto_combined_case1_thresh2, des_pareto_combined_case1_thresh2] = obtain_combined_pareto_data_case(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_feas_bools, constrad_read_case1, repeat3x3_case1, nfe_thresh2, sidenum, pop_size, num_runs); 

% Case 2
[f_pen_pareto_combined_case2_thresh2, f_true_pareto_combined_case2_thresh2, f_true_feas_pareto_combined_case2_thresh2, des_pareto_combined_case2_thresh2] = obtain_combined_pareto_data_case(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_feas_bools, constrad_read_case2, repeat3x3_case2, nfe_thresh2, sidenum, pop_size, num_runs); 

% NFE threshold = 1500
nfe_thresh3 = 1500;
% Case 1
[f_pen_pareto_combined_case1_thresh3, f_true_pareto_combined_case1_thresh3, f_true_feas_pareto_combined_case1_thresh3, des_pareto_combined_case1_thresh3] = obtain_combined_pareto_data_case(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_feas_bools, constrad_read_case1, repeat3x3_case1, nfe_thresh3, sidenum, pop_size, num_runs); 

% Case 2
[f_pen_pareto_combined_case2_thresh3, f_true_pareto_combined_case2_thresh3, f_true_feas_pareto_combined_case2_thresh3, des_pareto_combined_case2_thresh3] = obtain_combined_pareto_data_case(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_feas_bools, constrad_read_case2, repeat3x3_case2, nfe_thresh3, sidenum, pop_size, num_runs); 

% NFE threshold = 3000
nfe_thresh4 = 3000;
% Case 1
[f_pen_pareto_combined_case1_thresh4, f_true_pareto_combined_case1_thresh4, f_true_feas_pareto_combined_case1_thresh4, des_pareto_combined_case1_thresh4] = obtain_combined_pareto_data_case(fibre_model_used, case1_partcoll_bools, case1_nodalprop_bools, case1_orient_bools, case1_feas_bools, constrad_read_case1, repeat3x3_case1, nfe_thresh4, sidenum, pop_size, num_runs); 

% Case 2
[f_pen_pareto_combined_case2_thresh4, f_true_pareto_combined_case2_thresh4, f_true_feas_pareto_combined_case2_thresh4, des_pareto_combined_case2_thresh4] = obtain_combined_pareto_data_case(fibre_model_used, case2_partcoll_bools, case2_nodalprop_bools, case2_orient_bools, case2_feas_bools, constrad_read_case2, repeat3x3_case2, nfe_thresh4, sidenum, pop_size, num_runs); 


% Plotting
figure
subplot(2,2,1)
scatter(-f_pen_pareto_combined_case1_thresh1(:,1), f_pen_pareto_combined_case1_thresh1(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(-f_pen_pareto_combined_case2_thresh1(:,1), f_pen_pareto_combined_case2_thresh1(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold off
xlabel('Penalized $C_{22}$','Interpreter','Latex')
ylabel('Penalized $v_f$','Interpreter','Latex')
legend('Without Heuristics','With Heuristics','Location','best')
title('250 NFE')

subplot(2,2,2)
scatter(-f_pen_pareto_combined_case1_thresh2(:,1), f_pen_pareto_combined_case1_thresh2(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(-f_pen_pareto_combined_case2_thresh2(:,1), f_pen_pareto_combined_case2_thresh2(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold off
xlabel('Penalized $C_{22}$','Interpreter','Latex')
ylabel('Penalized $v_f$','Interpreter','Latex')
%legend('Without Heuristics','With Heuristics','Location','best')
title('500 NFE')

subplot(2,2,3)
scatter(-f_pen_pareto_combined_case1_thresh3(:,1), f_pen_pareto_combined_case1_thresh3(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(-f_pen_pareto_combined_case2_thresh3(:,1), f_pen_pareto_combined_case2_thresh3(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold off
xlabel('Penalized $C_{22}$','Interpreter','Latex')
ylabel('Penalized $v_f$','Interpreter','Latex')
%legend('Without Heuristics','With Heuristics','Location','best')
title('1500 NFE')

subplot(2,2,4)
scatter(-f_pen_pareto_combined_case1_thresh4(:,1), f_pen_pareto_combined_case1_thresh4(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(-f_pen_pareto_combined_case2_thresh4(:,1), f_pen_pareto_combined_case2_thresh4(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold off
xlabel('Penalized $C_{22}$','Interpreter','Latex')
ylabel('Penalized $v_f$','Interpreter','Latex')
%legend('Without Heuristics','With Heuristics','Location','best')
title('3000 NFE')

% Plotting feasible designs in true objectives space
figure
subplot(2,2,1)
scatter(f_true_feas_pareto_combined_case1_thresh1(:,1), f_true_feas_pareto_combined_case1_thresh1(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case2_thresh1(:,1), f_true_feas_pareto_combined_case2_thresh1(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold off
xlabel('$C_{22}$','Interpreter','Latex')
ylabel('$v_f$','Interpreter','Latex')
legend('Without Heuristics','With Heuristics','Location','best')
title('500 NFE')

subplot(2,2,2)
scatter(f_true_feas_pareto_combined_case1_thresh2(:,1), f_true_feas_pareto_combined_case1_thresh2(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case2_thresh2(:,1), f_true_feas_pareto_combined_case2_thresh2(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold off
xlabel('$C_{22}$','Interpreter','Latex')
ylabel('$v_f$','Interpreter','Latex')
%legend('Without Heuristics','With Heuristics','Location','best')
title('1000 NFE')

subplot(2,2,3)
scatter(f_true_feas_pareto_combined_case1_thresh3(:,1), f_true_feas_pareto_combined_case1_thresh3(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case2_thresh3(:,1), f_true_feas_pareto_combined_case2_thresh3(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold off
xlabel('$C_{22}$','Interpreter','Latex')
ylabel('$v_f$','Interpreter','Latex')
%legend('Without Heuristics','With Heuristics','Location','best')
title('3000 NFE')

subplot(2,2,4)
scatter(f_true_feas_pareto_combined_case1_thresh4(:,1), f_true_feas_pareto_combined_case1_thresh4(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'black') 
hold on
scatter(f_true_feas_pareto_combined_case2_thresh4(:,1), f_true_feas_pareto_combined_case2_thresh4(:,2), 'Marker', '+', 'MarkerEdgeColor', 'black') 
hold off
xlabel('$C_{22}$','Interpreter','Latex')
ylabel('$v_f$','Interpreter','Latex')
%legend('Without Heuristics','With Heuristics','Location','best')
title('10000 NFE')

%% Functions 
function [utopia_penobjs, utopia_trueobjs] = find_utopia_points(f_pen_pareto_case_thresh, f_true_pareto_case_thresh)
    
    
end

function [objs_pen_pareto_combined, objs_true_pareto_combined, objs_true_feas_pareto_combined, designs_pareto_combined] = obtain_combined_pareto_data_case(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_feas_bools, constrad_prob_read, repeat_case3x3, nfe_threshold, sidenodenum, n_pop, n_runs)
    
    n_members_total = nchoosek(sidenodenum^2,2);    
    data_array_struct = struct;
    designs_array_struct = struct;
    for i = 1:n_runs
        [data_array_case, designs_array_case] = read_csv_data_tillnfe(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, case_feas_bools, constrad_prob_read, repeat_case3x3, nfe_threshold, n_members_total, i-1);
        current_field = strcat('run_',num2str(i));
        data_array_struct.(current_field) = data_array_case;
        designs_array_struct.(current_field) = designs_array_case;
    end
%     if constrad_prob_read
%         for i = 1:n_runs
%             [data_array, designs_array = read_csv_data_tillnfe(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, constrad_prob_read, repeat_case3x3, nfe_threshold, n_members_total, i-1);
%         end
%     else
%         for i = 1:n_runs
%             [data_array, designs_array] = read_csv_data_tillnfe(fib_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, constrad_prob_read, repeat_case3x3, nfe_threshold, n_members_total, i-1);
%         end
%     end
    [pen_obj1_combined, pen_obj2_combined, true_obj1_combined, true_obj2_combined, feasvals_combined, des_combined] = create_combined_arrays(data_array_struct, designs_array_struct, constrad_prob_read, n_members_total, n_runs);
    f_pen_combined = [pen_obj1_combined, pen_obj2_combined];
    pareto_bool = paretofront(f_pen_combined);
    objs_pen_pareto_combined = f_pen_combined(pareto_bool==1,:);
    objs_true_pareto_combined = [true_obj1_combined(pareto_bool==1), true_obj2_combined(pareto_bool==1)];
    feasvals_pareto_combined = feasvals_combined(pareto_bool==1,:);
    designs_pareto_combined = des_combined(pareto_bool==1,:);
    objs_true_feas_pareto_combined = objs_true_pareto_combined(feasvals_pareto_combined==1,:);
end

function [data_array_req, design_array_req] = read_csv_data_tillnfe(fib_stiff, partcoll_bools, nodalprop_bools, orient_bools, feas_bools, constrad_read, repeat3x3_bool, nfe_to_reach, n_total_members, run_num)
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
            filename2 = strcat(filename2,"_prob2_fibre_fullpop.csv");
        else 
            filename2 = strcat(filename2,"_fibre_varrad_fullpop.csv");
        end
    else
        filepath3 = "Truss Stiffness\\";
        if constrad_read
            filename2 = strcat(filename2,"_prob2_truss_fullpop.csv");
        else
            filename2 = strcat(filename2,"_truss_varrad_fullpop.csv");
        end
    end
    
    %%%% read appropriate file 
    full_filepath = strcat(filepath,filepath_repeatcase,filepath_constrad,filepath2,filepath_moea,filepath3,filename,num2str(run_num),filename2);
    
    if constrad_read
        format_string = '%s';
        for j = 1:11
            format_string = strcat(format_string,'%f');
        end
        data_table = readtable(full_filepath,'Format',format_string,'HeaderLines',1);
    else
        format_string = '';
        for j = 1:(n_total_members+11)
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
    csv_data = data_table(:,end-10:end);
    
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

function [pen_obj1_combined, pen_obj2_combined, true_obj1_combined, true_obj2_combined, feas_combined, designs_combined] = create_combined_arrays(data_struct, designs_struct, read_constrad, n_total_members, n_runs)
    n_total = 0;
    for i = 1:n_runs
        current_field = strcat('run_',num2str(i));
        current_data_array = data_struct.(current_field);
        data_array_nonans_bool = any(isnan(current_data_array),2);
        data_array_nonans = current_data_array(~data_array_nonans_bool,:);
        n_total = n_total + size(data_array_nonans(:,1),1);
    end
    pen_obj1_combined = zeros(n_total,1);
    pen_obj2_combined = zeros(n_total,1);
    true_obj1_combined = zeros(n_total,1);
    true_obj2_combined = zeros(n_total,1);
    feas_combined = zeros(n_total,1);
    if read_constrad
        designs_combined = strings(n_total,1);
    else
        designs_combined = zeros(n_total,n_total_members);
    end
    index = 1;
    for i = 1:n_runs
        current_field = strcat('run_',num2str(i));
        current_data_array = data_struct.(current_field);
        data_array_nonans_bool = any(isnan(current_data_array),2);
        data_array_nonans = current_data_array(~data_array_nonans_bool,:);
        n_current = size(data_array_nonans(:,1),1);
        pen_obj1_combined(index:index+n_current-1,1) = data_array_nonans(:,2);
        pen_obj2_combined(index:index+n_current-1,1) = data_array_nonans(:,3);
        true_obj1_combined(index:index+n_current-1,1) = data_array_nonans(:,4);
        true_obj2_combined(index:index+n_current-1,1) = data_array_nonans(:,5);
        feas_combined(index:index+n_current-1,1) = data_array_nonans(:,6);
        if read_constrad
            current_designs_array = designs_struct.(current_field);
            design_array_nonans = current_designs_array(~data_array_nonans_bool,:);
            designs_combined(index:index+n_current-1,1) = design_array_nonans(:,1);
        else
            current_designs_array = designs_struct.(current_field);
            design_array_nonans = current_designs_array(~data_array_nonans_bool,:);
            designs_combined(index:index+n_current-1,:) = design_array_nonans;
        end
             
        index = index + n_current;
    end
end
