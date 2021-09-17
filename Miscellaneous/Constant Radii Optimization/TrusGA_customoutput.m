(%%
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
n_var = size(CA_all,1) - 4; % repeatable designs only 

%%% Running with population type "bitstring"
fitnessfcn = @(x)multiobjective_bitstring(x,CA_all,NC,A,E,sel,r,target_c_ratio); % objective function

% With custom output function
mutation_rate = 0.25;
ga_options = optimoptions('gamultiobj', 'PopulationType', 'bitstring', 'CrossoverFraction', 0.85, 'MutationFcn', {@mutationuniform, mutation_rate}, 'PopulationSize', 1000, 'FunctionTolerance', 1e-8, 'OutputFcn',@outputfcn_savegen2);
[x_pop, f_pop] = gamultiobj(fitnessfcn, n_var, [], [], [], [], [], [], ga_options); 

%% Find number of viable designs in each generation
n_gen = size(gagenarray,3);
n_viable = zeros(n_gen,1);
for i = 1:n_gen
    gen_array = gagenarray(:,:,i);
    viable_pos = find(~gen_array(:,2),1,'first');
    if isempty(viable_pos)
        viable_pos = 0;
    end
    n_viable(i,1) = viable_pos;
end

%% Plot designs in every 15 generations
f_vals = gagenarray(:,1:2,:);
pareto_bool_gens = plot_pareto_seak_customoutput(f_vals, n_viable);
