%% Truss MOEA
clear all
close all
clc

%% Defining Problem Constants
E = 10000; % Young's Modulus for polymeric material (example: 10000 Pa)
sel = 0.05; % Unit square side length (NOT individual truss length) (example: 5 cm)
r = 50*(10^-6); % Radius for cross-sectional area of (assumed circular) truss members (example: 50 micrometers)
A = pi*(r^2); % Cross-sectional area of truss member
% Nodal Coordinate Vector (Standard 3x3, 2D Grid) below (each row represents a node, first column is x-coordinates, second column is y-coordinates):
NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1]; 
CA_all = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 2,3; 2,4; 2,5; 2,6; 2,7; 2,8; 2,9; 3,4; 3,5; 3,6; 3,7; 3,8; 3,9; 4,5; 4,6; 4,7; 4,8; 4,9; 5,6; 5,7; 5,8; 5,9; 6,7; 6,8; 6,9; 7,8; 7,9; 8,9];
target_c_ratio = 1; % Ratio of C22/C11 to target
global CA_infeas;
CA_infeas = {};
global inf_con;
inf_con = 0;

param.CA_all = CA_all;
param.NC = NC;
param.A = A;
param.E = E;
param.sel = sel;
param.r = r;
param.c_ratio = target_c_ratio;
param.CA_infeasible = CA_infeas;
param.inf_count = inf_con;

%% Problem class parameters
clear eMOGA
eMOGA.objfun = 'multiobjective_evMOGA'; % m-function name for objectives computation
eMOGA.searchspaceUB = ones(1,size(CA_all,1)); % Search space upper bound
eMOGA.searchspaceLB = zeros(1,size(CA_all,1)); % Search space lower bound
eMOGA.objfun_dim = 2; % Objective space dimension
eMOGA.param = param; % Additional parameters for objective function
eMOGA.Nind_P = 500; % Population Size
eMOGA.Generations = 200; % Number of generations

%% Algorithm execution
[pfront,pset,eMOGA]=evMOGA(eMOGA);

%% Compute and plot feasibility and stability scores for final population
feas_scores = zeros(size(pset,1),1);
stab_scores = zeros(size(pset,1),1);
sidenum = 3;
for i = 1:size(pset,1)
    x_curr = pset(i,:);
    
    % Converting design vector to CA matrix
    %x_bin = x_curr>0.5;
    %CA_des = CA_all(x_bin~=0,:);
    
    CA_des = CA_all(x_curr~=0,:);
    
    %x_des = [0, 1, 0, 0, x_curr(1:2), 1, x_curr(3:4), 0, x_curr(5:18), ...
    %x_curr(19:end), 0, 1, 0]; % for forcing truss elements in 1,3; 1,7 and
    % 7;9
    %CA_des = CA_all(x_des~=0,:);
    
    %x_des = [x_curr(1), 0, x_curr(2:3), 0, 0, 0, 0, x_curr(4:7), 0, 0, 0, 0, ... 
        %x_curr(8:9), 0, 0, 0, x_curr(10), 0, x_curr(11:12), 0, x_curr(13:16), ...
        %0, x_curr(17:19), 0, x_curr(20)]; % only adjacent node connections allowed
    %CA_des = CA_all(x_des~=0,:);
    
    % Computing feasibility and stability scores
    feas_scores(i) = feasibility_checker_nonbinary(NC,CA_des);
    stab_scores(i) = stabilityTester_2D_updated(sidenum,CA_des,NC);
end

% Plotting
figure()
plot(feas_scores, stab_scores, 'b*')
xlabel('Feasibility scores')
ylabel('Stability scores')
title('Constraints comparison - Method 5')

%% Plotting Pareto Front (SEAK code)
%%% Finding only feasible designs
feasibility = false(size(pset,1),1);
feasibility_nonbinary = zeros(size(pset,1),1);
for i = 1:size(pset,1) 
    [feasibility_nonbinary(i), feasibility(i)] = feasibility_checker_boolean(pset(i,:), NC, CA_all); 
end
x_feas = pset(feasibility,:);
f_feas = pfront(feasibility,:);
%plot_pareto_seak(f_feas)

% Plotting the true feasible objectives
f_true_feas = zeros(size(x_feas,1),2);
lambda = 100;
for i = 1:size(x_feas,1)
    f_true_feas(i,:) = f_feas(i,:) + lambda*2*ones(1,2);  
end
plot_pareto_seak(f_true_feas,2)

%%% Plotting all designs
%plot_pareto_seak(pfront)

%% Convert design vector to binary array (not needed for bitstring design vector)
x_bin_feas = false(size(x_feas,1), size(CA_all,1));
for i = 1:size(x_feas,1)
    x_feas_i = x_feas(i,:);
    x_bin_feas(i,:) = x_feas_i>0.5; 
end

%% Visualize trusses
x_feas_unique = unique(x_feas,'rows');
for i = 1:size(x_feas_unique,1)
    visualize_truss_fromx_3x3(NC, CA_all, x_feas_unique(i,:))
end
