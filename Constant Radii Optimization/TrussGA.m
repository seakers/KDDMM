%%
clear all;
close all;
clc;

%% Defining problem constants
E = 10000; % Young's Modulus for polymeric material (example: 10000 Pa)
sel = 0.05; % Unit square side length (NOT individual truss length) (example: 5 cm)
r = 50*(10^-6); % Radius for cross-sectional area of (assumed circular) truss members (example: 50 micrometers)
A = pi*(r^2); % Cross-sectional area of truss member
% Nodal Coordinate Vector (Standard 3x3, 2D Grid) below (each row represents a node, first column is x-coordinates, second column is y-coordinates):
NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1]; 
CA_all = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 2,3; 2,4; 2,5; 2,6; 2,7; 2,8; 2,9; 3,4; 3,5; 3,6; 3,7; 3,8; 3,9; 4,5; 4,6; 4,7; 4,8; 4,9; 5,6; 5,7; 5,8; 5,9; 6,7; 6,8; 6,9; 7,8; 7,9; 8,9];
target_c_ratio = 1; % Ratio of C22/C11 to target

%% Genetic Algorithm
%n_var = size(CA_all,1); % size of design vector
%n_var = size(CA_all,1) - 9; % for forcing truss elements in 1,3; 1,7 and
%7;9
%n_var = size(CA_all,1) - 16; % only adjacent node connections allowed
n_var = size(CA_all,1) - 4; % repeatable designs only 
%global CA_infeas;
%CA_infeas = {};
%global inf_con;
%inf_con = 0;
%global lambda;
%lambda = 0.5;
%global iter;
%iter = 0;

%%% Running with population type "double vector" and constraints 
%lb = zeros(1,size(CA_all,1)); % lower limit of design vector values
%ub = 9*ones(1,size(CA_all,1)); % upper limit of design vector values
%constraint_function = @(x)constraint_checker_nonbinary(x,NC,CA_all,inf_con,CA_infeas); % non-linear constraint function
% Use custom creation function
%ga_options = optimoptions('gamultiobj', 'CreationFcn', @customCreateFunction, 'PopulationType', 'doubleVector', 'PopulationSize', 500, 'FunctionTolerance', 1e-8); %'PlotFcn', 'gaplotpareto');

%ga_options = optimoptions('gamultiobj', 'PopulationType', 'doubleVector', 'PopulationSize', 200, 'FunctionTolerance', 1e-8); %'PlotFcn', 'gaplotpareto');
%fitnessfcn = @(x)multiobjective(x,CA_all,NC,A,E,sel,r,target_c_ratio,CA_infeas,inf_con); % objective function
%[x_pop, f_pop] = gamultiobj(fitnessfcn, n_var, [], [], [], [], lb, ub, constraint_function, ga_options);

%IntCon = linspace(1,n_var,n_var);
%[x_pop, f_pop] = gamultiobj(fitnessfcn, n_var, [], [], [], [], lb, ub, constraint_function, IntCon, ga_options); 

%%% Running with population type "bitstring"
fitnessfcn = @(x)multiobjective_bitstring(x,CA_all,NC,A,E,sel,r,target_c_ratio); % objective function

% Use custom creation function
%ga_options = optimoptions('gamultiobj', 'CreationFcn', @customCreateFunction, 'PopulationType', 'bitstring', 'PopulationSize', 500, 'FunctionTolerance', 1e-8);

mutation_rate = 0.25;
ga_options = optimoptions('gamultiobj', 'PopulationType', 'bitstring', 'CrossoverFraction', 0.85, 'MutationFcn', {@mutationuniform, mutation_rate}, 'PopulationSize', 1000, 'FunctionTolerance', 1e-8);
[x_pop, f_pop] = gamultiobj(fitnessfcn, n_var, [], [], [], [], [], [], ga_options); 

%% Compute and plot feasibility and stability scores for final population
feas_scores = zeros(size(x_pop,1),1);
stab_scores = zeros(size(x_pop,1),1);
sidenum = 3;
for i = 1:size(x_pop,1)
    x_curr = x_pop(i,:);
    
    x_des = [x_curr(1), x_curr(2), x_curr(3), x_curr(4), x_curr(5), x_curr(6), ...
        x_curr(7), x_curr(8), x_curr(9), x_curr(10), x_curr(11), x_curr(12), ...
        x_curr(13), x_curr(14), x_curr(15), x_curr(16), x_curr(17), x_curr(3), ...
        x_curr(18), x_curr(19), x_curr(20), x_curr(21), x_curr(22), x_curr(23), ...
        x_curr(24), x_curr(25), x_curr(26), x_curr(27), x_curr(28), x_curr(29), ...
        x_curr(30), x_curr(31), x_curr(23), x_curr(1), x_curr(32), x_curr(9)]; 
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
title('Constraints comparison - Method 6')

%% Write infeasible designs to text file
%fileId = fopen('CA_infeasible.txt','w');
%fprintf(fileId, '%6s %12s\r\n', 'Sr.No.', 'Infeasible CA')
%for i = 1:inf_con
    %CA_inf_current = CA_infeas{inf_con};
    %fprintf(fileId, '%d\t', i)
    %for j = 1:size(CA_inf_current,1)
        %row = CA_inf_current(j,:);
        %for k = 1:size(CA_inf_current,2)
            %if k==size(CA_inf_current,2)
               %fprintf(fileId, '%d\n', row(k)) 
            %else
               %fprintf(fileId, '%d', row(k))
            %end
        %end
    %end
    %fprintf(fileId, '\n')
%end
%fclose(fileId);

%% Plotting Pareto Front (SEAK code)
%%% For 2 objectives
% Plot only feasible objective values

%feasibility = false(size(x_pop,1),1);
%stability = false(size(x_pop,1),1);
%feasibility_nonbinary = zeros(size(x_pop,1),1);
%stability_nonbinary = zeros(size(x_pop,1),1);
%sidenum = 3;
%for i = 1:size(x_pop,1) 
    %x_vec = x_pop(i,:);
    %x_des = [x_vec(1), x_vec(2), x_vec(3), x_vec(4), x_vec(5), x_vec(6), ...
        %x_vec(7), x_vec(8), x_vec(9), x_vec(10), x_vec(11), x_vec(12), ...
        %x_vec(13), x_vec(14), x_vec(15), x_vec(16), x_vec(17), x_vec(3), ...
        %x_vec(18), x_vec(19), x_vec(20), x_vec(21), x_vec(22), x_vec(23), ...
        %x_vec(24), x_vec(25), x_vec(26), x_vec(27), x_vec(28), x_vec(29), ...
        %x_vec(30), x_vec(31), x_vec(23), x_vec(1), x_vec(32), x_vec(9)]; 
        % repeatable designs
    %CA_binary = CA_all(x_des~=0,:);
    %[feasibility_nonbinary(i), feasibility(i)] = feasibility_checker_boolean(x_des, NC, CA_all); 
    %%stability(i) = stabilityTester_2D_boolean(sidenum,x_des,CA_all,NC);
    %stability_nonbinary(i) = stabilityTester_2D_updated(sidenum, CA_binary, NC); 
%end

x_feas = x_pop(((feas_scores==1)&(stab_scores>=0.8)),:);
f_feas = f_pop(((feas_scores==1)&(stab_scores>=0.8)),:);

%x_feas = x_pop(feasibility&(stability>=0.8),:);
%f_feas = f_pop(feasibility&(stability>=0.8),:);
%plot_pareto_seak(f_feas,2)

% Plotting the true feasible objectives
f_true_feas = zeros(size(x_feas,1),2);
stab_feas = stab_scores(((feas_scores==1)&(stab_scores>=0.8)),:);
%stab_feas = false(size(x_feas,1),1);
%alpha = 1e-4;
for i = 1:size(x_feas,1)
    %x_vec = x_feas(i,:);
    %x_des = [x_vec(1), x_vec(2), x_vec(3), x_vec(4), x_vec(5), x_vec(6), ...
        %x_vec(7), x_vec(8), x_vec(9), x_vec(10), x_vec(11), x_vec(12), ...
        %x_vec(13), x_vec(14), x_vec(15), x_vec(16), x_vec(17), x_vec(3), ...
        %x_vec(18), x_vec(19), x_vec(20), x_vec(21), x_vec(22), x_vec(23), ...
        %x_vec(24), x_vec(25), x_vec(26), x_vec(27), x_vec(28), x_vec(29), ...
        %x_vec(30), x_vec(31), x_vec(23), x_vec(1), x_vec(32), x_vec(9)]; 
        % repeatable designs
    %CA_feas = CA_all(x_des~=0,:);
    %stab_feas = stabilityTester_2D_updated(sidenum, CA_feas, NC);
    
    stab_feas_des = stab_feas(i);

    %lambda = 0.5*(1 - (alpha*i));
    %lambda1 = 0.5;
    %lambda2 = 6000;
    
    f_true_feas(i,:) = [15*(f_feas(i,1) + log(abs(stab_feas_des))/2), 6000*(f_feas(i,2) + log(abs(stab_feas_des))/2)];
    %f_true_feas(i,:) = [15*(f_feas(i,1) + (tanh(1.1)+tanh(stab_feas_des+0.1))/2), 6000*(f_feas_des(i,2) + (tanh(1.1)+tanh(stab_feas_des+0.1))/2)];
    %f_true_feas(i,:) = [15*(f_feas(i,1) + (1+stab_feas_des)/2), 6000*(f_feas(i,2) + (1+stab_feas_des)/2)];  
end
figure
pareto_bool = plot_pareto_seak(f_true_feas,2);

% Plot all objective values
%figure
%plot_pareto_seak(f_pop,2)

%%% For 3 objectives
%figure
%plot_pareto_seak(f_pop,3)

%% Convert feasible and stable design vector to binary array (not needed for binary designs)
%x_bin_feas = false(size(x_feas,1), size(CA_all,1));
%for i = 1:size(x_feas,1)
    %x_feas_i = x_feas(i,:);
    %x_bin_feas(i,:) = x_feas_i>0.5; 
%end
%x_feas = x_bin_feas;

%% Visualize pareto designs 
x_pareto = x_feas(pareto_bool==1,:);
x_pareto_unique = unique(x_pareto,'rows');
for i = 1:size(x_pareto_unique,1)
    visualize_truss_fromx_3x3(NC, CA_all, x_pareto_unique(i,:))
end

%% check neighbouring designs
%scale_mat = [15,0;0,6000];
n_des_var = size(CA_all,1);
x_pareto_unique = unique(x_feas,'rows');
f_des = zeros(1,2);
f_pareto_unique = zeros(size(x_pareto_unique,1),2);
for i = 1:size(x_pareto_unique,1)
    x_vec = x_pareto_unique(i,:);
    x_des = [x_vec(1), x_vec(2), x_vec(3), x_vec(4), x_vec(5), x_vec(6), ...
        x_vec(7), x_vec(8), x_vec(9), x_vec(10), x_vec(11), x_vec(12), ...
        x_vec(13), x_vec(14), x_vec(15), x_vec(16), x_vec(17), x_vec(3), ...
        x_vec(18), x_vec(19), x_vec(20), x_vec(21), x_vec(22), x_vec(23), ...
        x_vec(24), x_vec(25), x_vec(26), x_vec(27), x_vec(28), x_vec(29), ...
        x_vec(30), x_vec(31), x_vec(23), x_vec(1), x_vec(32), x_vec(9)]; 
        % repeatable designs
    CA_feas_unique = CA_all(x_des~=0,:);
    stab_feas_unique = stabilityTester_2D_updated(sidenum, CA_feas_unique, NC);
    f_des = multiobjective_bitstring(x_vec,CA_all,NC,A,E,sel,r,target_c_ratio);
    
    f_pareto_unique(i,:) = [15*(f_des(1) + log(abs(stab_feas_unique))/2), 6000*(f_des(2) + log(abs(stab_feas_unique))/2)];
    %f_pareto_unique(i,:) = [15*(f_des(1) + (tanh(1.1)+tanh(stab_feas_unique+0.1))/2), 6000*(f_des(2) + (tanh(1.1)+tanh(stab_feas_unique+0.1))/2)];
    %f_pareto_unique(i,:) = [15*(f_des(1) + (1+stab_feas_unique)/2), 6000*(f_des(2) + (1+stab_feas_unique)/2)];   
end

feas_and_stab = false(size(x_pareto_unique(1,:),1),size(x_pareto_unique(1,:),2));
f_des_mod = zeros(1,2);
n_feas_and_stab = zeros(size(x_pareto_unique(1,:),1),1);
x_feas_and_stab = false((size(x_pareto_unique,1)*n_des_var),size(x_pareto_unique,2));
f_feas_and_stab = zeros((size(x_pareto_unique,1)*n_des_var),2);
constraint_scores = zeros((size(x_pareto_unique,1)*n_des_var),2);
count = 0;

%%% Check the neighbourhood of pareto set and store the feasible and stable
%%% neighbourhood designs and evaluations
for i = 1:size(x_pareto_unique,1)
    %des = x_feas_unique(i,:);
    
    x_vec = x_pareto_unique(i,:);
    des = [x_vec(1), x_vec(2), x_vec(3), x_vec(4), x_vec(5), x_vec(6), ...
        x_vec(7), x_vec(8), x_vec(9), x_vec(10), x_vec(11), x_vec(12), ...
        x_vec(13), x_vec(14), x_vec(15), x_vec(16), x_vec(17), x_vec(3), ...
        x_vec(18), x_vec(19), x_vec(20), x_vec(21), x_vec(22), x_vec(23), ...
        x_vec(24), x_vec(25), x_vec(26), x_vec(27), x_vec(28), x_vec(29), ...
        x_vec(30), x_vec(31), x_vec(23), x_vec(1), x_vec(32), x_vec(9)]; 
    % repeatable designs  
    
    %constraint_scores_des = zeros(n_des_var,2);
    for j = 1:size(des,2)
        des_mod = des;
        des_mod(j) = ~des(j);
        CA_des_mod = CA_all(des_mod~=0,:);
        
        feas_des_mod = feasibility_checker_nonbinary(NC,CA_des_mod);
        stab_des_mod = stabilityTester_2D_updated(sidenum,CA_des_mod,NC);
        
        %constraint_scores_des(j,:) = [feas_des_mod, stab_des_mod];
        
        if ((feas_des_mod==1) && (stab_des_mod>=0.8))
            feas_and_stab(i,j) = true;
            x_des_mod = [des_mod(1:17),des_mod(19:32),des_mod(35)]; 
            x_feas_and_stab(count+1,:) = x_des_mod;
            f_des_mod = multiobjective_bitstring(x_des_mod,CA_all,NC,A,E,sel,r,target_c_ratio);
            
            f_feas_and_stab(count+1,:) = [15*(f_des_mod(1) + log(abs(stab_des_mod))/2), 6000*(f_des_mod(2) + log(abs(stab_des_mod))/2)];
            %f_feas_and_stab(count+1,:) = [15*(f_des_mod(1) + (tanh(1.1)+tanh(stab_des_mod+0.1))/2), 6000*(f_des_mod(2) + (tanh(1.1)+tanh(stab_des_mod+0.1))/2)];
            %f_feas_and_stab(count+1,:) = [15*(f_des_mod(1) + (1+stab_des_mod)/2), 6000*(f_des_mod(2) + (1+stab_des_mod)/2)]; 
            
            constraint_scores(count+1,:) = [feas_des_mod, stab_des_mod]; 
            
            count = count + 1;
        end
    end
    n_feas_and_stab(i,1) = nnz(feas_and_stab(i,:));
    %constraint_scores(((i-1)*n_des_var)+1:n_des_var*i,:) = constraint_scores_des; 
end
f_combined = vertcat(f_pareto_unique, f_feas_and_stab(1:count,:));
x_combined = vertcat(x_pareto_unique, x_feas_and_stab(1:count,:));
n_total_viable = sum(n_feas_and_stab);

%%% Plot the pareto set with the feasible and stable neighbourhood
%%% evaluations
figure
par_bool = plot_pareto_seak(f_combined,2);
x_pareto_combined = x_combined(par_bool==1,:);

%pareto_bool1 = paretofront(f_feas_unique);
%f_pareto1 = f_feas_unique(pareto_bool1==1,:);
%f_pareto_true1 = [f_pareto1(:,1), -f_pareto1(:,2)]; % objective 1 is to be
% minimized while objective 2 is to be maximized 
    
%f_additional = f_feas_and_stab(1:count,:);
%pareto_bool2 = paretofront(f_additional);
%f_pareto2 = f_additional(pareto_bool2==1,:);
%f_pareto_true2 = [f_pareto2(:,1), -f_pareto2(:,2)]; % objective 1 is to be
% minimized while objective 2 is to be maximized 

%figure
%plot(f_pareto_true1(:,1), f_pareto_true1(:,2),'*b')
%hold on 
%plot(f_pareto_true2(:,1), f_pareto_true2(:,2),'*r')
%hold off
%xlabel('$\left|\frac{C_{22}}{C_{11}} - c_{target}\right|$','Interpreter','latex')
%ylabel('$\frac{C_{22}}{v_f}$','Interpreter','latex')
