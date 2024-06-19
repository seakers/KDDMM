%% Check ease of satisfaction of different heuristics
clear
close all
clc

%% Parameter definitions (printable designs)
E = 1e6; % Young's Modulus for polymeric material (example: 1 MPa)
sel = 0.01; % Unit square side length (NOT individual truss length) (in m)
r = 25*(10^-5); % Radius for cross-sectional area of (assumed circular) truss members (in m)
A = pi*(r^2); % Cross-sectional area of truss member
NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1]; 
CA_all = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 2,3; 2,4; 2,5; 2,6; 2,7; 2,8; 2,9; 3,4; 3,5; 3,6; 3,7; 3,8; 3,9; 4,5; 4,6; 4,7; 4,8; 4,9; 5,6; 5,7; 5,8; 5,9; 6,7; 6,8; 6,9; 7,8; 7,9; 8,9];
c_ratio = 1;
sidenum = 3;

%% Check ease of satisfaction for Feasibility
% Generate random architectures 
n_des = 500; % number of architectures to generate
n_runs = 10; % number of runs to generate "n_des" architectures
bool_des_all = zeros(36,n_des,n_runs);
feas_all = zeros(n_des,n_runs);
n_full_feas_runs = zeros(n_runs,1);
n_good_feas_runs = zeros(n_runs,1);
n_nz_feas_runs = zeros(n_runs,1);
for i = 1:n_runs
    feas_run = zeros(n_des,1);
    bool_des_run = zeros(36,n_des);
    for j = 1:n_des
        bool_des = randi([0,1],32,1);
        complete_bool_des = [bool_des(1:17); bool_des(3); bool_des(18:31); bool_des(23); bool_des(1); bool_des(32); bool_des(9)];
        CA_des = CA_all(complete_bool_des~=0,:);
        bool_des_run(:,j) = complete_bool_des;
        feas_run(j) = feasibility_checker_nonbinary_V2(NC,CA_des);
    end
    feas_all(:,i) = feas_run;
    bool_des_all(:,:,i) = bool_des_run;
    n_full_feas_runs(i,1) = length(feas_run(feas_run==1));
    n_good_feas_runs(i,1) = length(feas_run(feas_run>=0.8));
    n_nz_feas_runs(i,1) = length(feas_run(feas_run>=0.1));
end
ratio_full_feas_runs = n_full_feas_runs./n_des;
mean_ratio_full_feas = mean(ratio_full_feas_runs)

ratio_good_feas_runs = n_good_feas_runs./n_des;
mean_ratio_good_feas = mean(ratio_good_feas_runs)

ratio_nz_feas_runs = n_nz_feas_runs./n_des;
mean_ratio_nz_feas = mean(ratio_nz_feas_runs)

%% Check ease of satisfaction for Stability Heuristic
% Generate random architectures 
n_des = 500; % number of architectures to generate
n_runs = 10; % number of runs to generate "n_des" architectures
bool_des_all = zeros(36,n_des,n_runs);
stab_all = zeros(n_des,n_runs);
n_full_stab_runs = zeros(n_runs,1);
n_good_stab_runs = zeros(n_runs,1);
n_nz_stab_runs = zeros(n_runs,1);
for i = 1:n_runs
    stab_run = zeros(n_des,1);
    bool_des_run = zeros(36,n_des);
    for j = 1:n_des
        bool_des = randi([0,1],32,1);
        complete_bool_des = [bool_des(1:17); bool_des(3); bool_des(18:31); bool_des(23); bool_des(1); bool_des(32); bool_des(9)];
        CA_des = CA_all(complete_bool_des~=0,:);
        bool_des_run(:,j) = complete_bool_des;
        stab_run(j) = stabilityTester_2D_V7(sidenum,CA_des,NC,sel);
    end
    stab_all(:,i) = stab_run;
    bool_des_all(:,:,i) = bool_des_run;
    n_full_stab_runs(i,1) = length(stab_run(stab_run==1));
    n_good_stab_runs(i,1) = length(stab_run(stab_run>=0.8));
    n_nz_stab_runs(i,1) = length(stab_run(stab_run>=0.1));
end
ratio_full_stab_runs = n_full_stab_runs./n_des;
mean_ratio_full_stab = mean(ratio_full_stab_runs)

ratio_good_stab_runs = n_good_stab_runs./n_des;
mean_ratio_good_stab = mean(ratio_good_stab_runs)

ratio_nz_stab_runs = n_nz_stab_runs./n_des;
mean_ratio_nz_stab = mean(ratio_nz_stab_runs)

%% Check ease of satisfaction for Orientation Heuristic
% Generate random architectures 
n_des = 500; % number of architectures to generate
n_runs = 10; % number of runs to generate "n_des" architectures
bool_des_all = zeros(36,n_des,n_runs);
orient_all = zeros(n_des,n_runs);
n_full_orient_runs = zeros(n_runs,1);
n_good_orient_runs = zeros(n_runs,1);
n_nz_orient_runs = zeros(n_runs,1);
for i = 1:n_runs
    orient_run = zeros(n_des,1);
    bool_des_run = zeros(36,n_des);
    for j = 1:n_des
        bool_des = randi([0,1],32,1);
        complete_bool_des = [bool_des(1:17); bool_des(3); bool_des(18:31); bool_des(23); bool_des(1); bool_des(32); bool_des(9)];
        CA_des = CA_all(complete_bool_des~=0,:);
        bool_des_run(:,j) = complete_bool_des;
        [orient_run(j), ~] = orientationHeuristic_V2(NC,CA_des,c_ratio);
    end
    orient_all(:,i) = orient_run;
    bool_des_all(:,:,i) = bool_des_run;
    n_full_orient_runs(i,1) = length(orient_run(orient_run==1));
    n_good_orient_runs(i,1) = length(orient_run(orient_run>=0.8));
    n_nz_orient_runs(i,1) = length(orient_run(orient_run>=0.1));
end
ratio_full_orient_runs = n_full_orient_runs./n_des;
mean_ratio_full_orient = mean(ratio_full_orient_runs)

ratio_good_orient_runs = n_good_orient_runs./n_des;
mean_ratio_good_orient = mean(ratio_good_orient_runs)

ratio_nz_orient_runs = n_nz_orient_runs./n_des;
mean_ratio_nz_orient = mean(ratio_nz_orient_runs)

%% Check ease of satisfaction for Printability Heuristic

% Compute mean satisfaction ratios for different "sel" values
sel_array = linspace(0.004,0.01,7);
mean_full_ratios = zeros(length(sel_array),1);
mean_good_ratios = zeros(length(sel_array),1);
mean_nz_ratios = zeros(length(sel_array),1);
for k = 1:length(sel_array)
    [mean_ratio_full_current,mean_ratio_good_current,mean_ratio_nz_current] = find_eos_print(sel_array(k),r,CA_all);
    mean_full_ratios(k,1) = mean_ratio_full_current;
    mean_good_ratios(k,1) = mean_ratio_good_current;
    mean_nz_ratios(k,1) = mean_ratio_nz_current;
end

% Plotting
figure
plot(sel_array',mean_full_ratios.*100,'*b');
xlabel('Side Element Length in m')
ylabel('Mean % of satisfaction (10 runs)')

figure
plot(sel_array',mean_good_ratios.*100,'*b');
xlabel('Side Element Length in m')
ylabel('Mean % of good designs (10 runs)')

figure
plot(sel_array',mean_nz_ratios.*100,'*b');
xlabel('Side Element Length in m')
ylabel('Mean % of nonzero printability designs (10 runs)')

function [mean_ratio_full_print, mean_ratio_good_print, mean_ratio_nz_print] = find_eos_print(sidel,rad,CA_full)
    % Define Nodal Connectivity array
    NC = sidel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1]; 
    
    % Generate random architectures 
    n_des = 500; % number of architectures to generate
    n_runs = 10; % number of runs to generate "n_des" architectures
    bool_des_all = zeros(36,n_des,n_runs);
    print_all = zeros(n_des,n_runs);
    n_full_print_runs = zeros(n_runs,1);
    n_good_print_runs = zeros(n_runs,1);
    n_nz_print_runs = zeros(n_runs,1);
    for i = 1:n_runs
        print_run = zeros(n_des,1);
        bool_des_run = zeros(36,n_des);
        for j = 1:n_des
            bool_des = randi([0,1],32,1);
            complete_bool_des = [bool_des(1:17); bool_des(3); bool_des(18:31); bool_des(23); bool_des(1); bool_des(32); bool_des(9)];
            CA_des = CA_full(complete_bool_des~=0,:);
            bool_des_run(:,j) = complete_bool_des;
            rvar = rad.*ones(1,size(CA_des,1));
            [print_run(j), ~] = printChecker(NC,CA_des,rvar);
        end
        print_all(:,i) = print_run;
        bool_des_all(:,:,i) = bool_des_run;
        n_full_print_runs(i,1) = length(print_run(print_run==1));
        n_good_print_runs(i,1) = length(print_run(print_run>=0.8));
        n_nz_print_runs(i,1) = length(print_run(print_run>=0.1));
    end
    ratio_full_print_runs = n_full_print_runs./n_des;
    mean_ratio_full_print = mean(ratio_full_print_runs);
    
    ratio_good_print_runs = n_good_print_runs./n_des;
    mean_ratio_good_print = mean(ratio_good_print_runs);
    
    ratio_nz_print_runs = n_nz_print_runs./n_des;
    mean_ratio_nz_print = mean(ratio_nz_print_runs);
    
    % Find % of non-zero printability designs
    %print_all_bool = print_all>=0.1;
    %nnz_print_all = nnz(print_all_bool);
    %percent_nnz_print_all = nnz_print_all*100/(n_des*n_runs);
end