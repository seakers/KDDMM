%% Plot comparison between AOS and Epsilon MOEA
clear
close all
clc

%% Read results csv file
%%%% set to: 
%       true - to read fibre stiffness evaluation run results
%       false - to read truss stiffness evaluation run results
fibre_stiffness1 = true;
fibre_stiffness2 = false;

%%%% set to: 
%       true - to read Epsilon MOEA run results
%       false - to read AOS MOEA run results
epsilon_moea1 = true;
epsilon_moea2 = false;

%%%% set to: 
%       0 - both feasibility and stability is not used in both operators
%       1 - only feasibility is used ine AddTruss operators
%       2 - both feasibility and stability are used in both operators
feas_and_stab = 0;

%%%% read appropriate files 
full_filepath_eps = extract_filename(fibre_stiffness2, epsilon_moea1, feas_and_stab);
data_table_eps = readtable(full_filepath_eps,'Format','%s%f%f%f%f','HeaderLines',1);

full_filepath_aos = extract_filename(fibre_stiffness2, epsilon_moea2, feas_and_stab);
data_table_aos = readtable(full_filepath_aos,'Format','%s%f%f%f%f','HeaderLines',1);

%%%% store retrieved data into different variables
%%%% csv_data includes: [Pen. Obj. 1, Pen.Obj. 2, Feasibility Score,
%%%% Stablity Score]
pop_size_eps =  size(data_table_eps,1);
csv_data_eps = zeros(pop_size_eps,4);
designs_eps = strings(pop_size_eps);
csv_data_eps = data_table_eps(:,2:end);
designs_eps = data_table_eps(:,1);

pop_size_aos =  size(data_table_aos,1);
csv_data_aos = zeros(pop_size_aos,4);
designs_aos = strings(pop_size_aos);
csv_data_aos = data_table_aos(:,2:end);
designs_aos = data_table_aos(:,1);

%% Compute true objectives
f_true_eps = compute_true_objectives(csv_data_eps, pop_size_eps, fibre_stiffness2);
f_true_aos = compute_true_objectives(csv_data_aos, pop_size_aos, fibre_stiffness2);

%% Plotting
plot_pareto_seak_compare(f_true_eps,f_true_aos);
