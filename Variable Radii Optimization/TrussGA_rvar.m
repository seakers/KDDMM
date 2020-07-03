%% Truss GA rvar
clear 
close all
clc

%% Defining problem parameters
E = 10000; % Young's Modulus for polymeric material (example: 10000 Pa)
sel = 0.05; % Unit square side length (NOT individual truss length) (example: 5 cm)
CA_all = [1,2; 1,4; 1,5; 2,3; 2,4; 2,5; 2,6; 3,5; 3,6; 4,5; 4,7; 4,8; 5,6; 5,7; 5,8; 5,9; 6,8; 6,9; 7,8; 8,9]; % only adjacent nodes
nucFac = 1;
sidenum = (2*nucFac) + 1; 
NC = generateNC(sel, sidenum);
target_c_ratio = 1; % Ratio of C22/C11 to target

%% Genetic algorithm
n_var = size(CA_all,1) - 4; % repeatable designs only
lb = zeros(1,n_var); % lower bound for radii values
ub = 1e-4.*ones(1,n_var); % upper bound for radii values

% With interior penalty functions
%fitnessfcn = @(x)multiobjective_rvar(x, nucFac, sel, E, CA_all, target_c_ratio);
%mutation_rate = 0.25;
%ga_options = optimoptions('gamultiobj', 'PopulationType', 'doubleVector', 'CrossoverFraction', 0.85, 'MutationFcn', @mutationadaptfeasible, 'PopulationSize', 1000, 'FunctionTolerance', 1e-8);
%[x_pop, f_pop] = gamultiobj(fitnessfcn, n_var, [], [], [], [], lb, ub, ga_options); 

% With nonlinear constraint function
fitnessfcn = @(x)multiobjective_rvar_nonlcon(x, nucFac, sel, E, CA_all, target_c_ratio);
constraintfcn = @(x)constraint_checker_rvar(x, nucFac, sel, CA_all); 
%mutation_rate = 0.25;
ga_options = optimoptions('gamultiobj', 'PopulationType', 'doubleVector', 'CrossoverFraction', 0.85, 'MutationFcn', @mutationadaptfeasible, 'PopulationSize', 1000, 'FunctionTolerance', 1e-8);
[x_pop, f_pop] = gamultiobj(fitnessfcn, n_var, [], [], [], [], lb, ub, constraintfcn, ga_options); 

%% Compute and plot feasibility and stability scores for final population
feas_scores = zeros(size(x_pop,1),1);
stab_scores = zeros(size(x_pop,1),1);

for i = 1:size(x_pop,1)
    x_curr = x_pop(i,:);
    
    x_des = [x_curr(1:8), x_curr(2), x_curr(9:16), x_curr(11), x_curr(1), x_curr(4)];
    % repeatable designs
    CA_des = CA_all(x_des~=0,:);
    
    % Computing feasibility and stability scores
    feas_scores(i) = feasibility_checker_nonbinary(NC,CA_des);
    stab_scores(i) = stabilityTester_2D_updated(sidenum,CA_des,NC);
end

% Plotting
figure()
plot(feas_scores, stab_scores, 'b*')
xlabel('Feasibility scores')
ylabel('Stability scores')
title('Constraints comparison - Variable Radii')

%% Plotting Pareto Front (SEAK code)
% Plot only feasible objective values
x_feas = x_pop(((feas_scores==1)&(stab_scores>=0.8)),:);
f_feas = f_pop(((feas_scores==1)&(stab_scores>=0.8)),:);

%x_feas = x_pop(feasibility&(stability>=0.8),:);
%f_feas = f_pop(feasibility&(stability>=0.8),:);
%plot_pareto_seak(f_feas,2)

% Plotting the true feasible objectives
f_true_feas = zeros(size(x_feas,1),2);
stab_feas = stab_scores(((feas_scores==1)&(stab_scores>=0.8)),:);
%stab_feas = false(size(x_feas,1),1);
for i = 1:size(x_feas,1)  
    stab_feas_des = stab_feas(i);
    
    f_true_feas(i,:) = [15*(f_feas(i,1) + log(abs(stab_feas_des))/2), 6000*(f_feas(i,2) + log(abs(stab_feas_des))/2)];
end
figure
pareto_bool = plot_pareto_seak(f_true_feas,2);
