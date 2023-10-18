%% Plot combined pareto front comparison
clear all
close all
clc

%% Read results from csv files
%%%% set to: 
%       true - to read fibre stiffness evaluation run results
%       false - to read truss stiffness evaluation run results
fibre_stiffness = true;  

%%%% set algo_choice to: 
%       eps_moea - to read Epsilon MOEA run results
%       aos - to read AOS MOEA run results
%       sc_dnf - to read DNF run results
%       sc_ach - to read ACH run results
%algo_choice = 'sc_ach';

%%%% set to: 
%       0 - both feasibility and stability is not used in both operators
%       1 - only feasibility is used ine AddTruss operators
%       2 - both feasibility and stability are used in both operators
feas_and_stab = 0;

%%%% set to:
%       true - to read biased initialization results
%       false - to read random initialization results
biased_init = true;

num_runs = 30; % change based on run 
pop_size = 100; % change based on run

%% Generate combined data arrays for each run case
% For epsilon MOEA, fibre stiffness, random initialization
data_array_eps_fib_rand = zeros(pop_size,4,num_runs);
designs_array_eps_fib_rand = strings(pop_size,num_runs);
for i = 1:num_runs
    [data_array_eps_fib_rand(:,:,i), designs_array_eps_fib_rand(:,i)] = read_csv_file(fibre_stiffness, 'eps_moea', feas_and_stab, ~biased_init, i-1);
end
[obj1_combined_eps_fib_rand, obj2_combined_eps_fib_rand, feas_combined_eps_fib_rand, stab_combined_eps_fib_rand, des_combined_eps_fib_rand] = create_combined_arrays(data_array_eps_fib_rand, designs_array_eps_fib_rand, num_runs);

% For epsilon MOEA, truss stiffness, random initialization
data_array_eps_truss_rand = zeros(pop_size,4,num_runs);
designs_array_eps_truss_rand = strings(pop_size,num_runs);
for i = 1:num_runs
    [data_array_eps_truss_rand(:,:,i), designs_array_eps_truss_rand(:,i)] = read_csv_file(~fibre_stiffness, 'eps_moea', feas_and_stab, ~biased_init, i-1);
end
[obj1_combined_eps_truss_rand, obj2_combined_eps_truss_rand, feas_combined_eps_truss_rand, stab_combined_eps_truss_rand, des_combined_eps_truss_rand] = create_combined_arrays(data_array_eps_truss_rand, designs_array_eps_truss_rand, num_runs);

% For epsilon MOEA, fibre stiffness, biased initialization
data_array_eps_fib_bias = zeros(pop_size,4,num_runs);
designs_array_eps_fib_bias = strings(pop_size,num_runs);
for i = 1:num_runs
    [data_array_eps_fib_bias(:,:,i), designs_array_eps_fib_bias(:,i)] = read_csv_file(fibre_stiffness, 'eps_moea', feas_and_stab, biased_init, i-1);
end
[obj1_combined_eps_fib_bias, obj2_combined_eps_fib_bias, feas_combined_eps_fib_bias, stab_combined_eps_fib_bias, des_combined_eps_fib_bias] = create_combined_arrays(data_array_eps_fib_bias, designs_array_eps_fib_bias, num_runs);

% For epsilon MOEA, truss stiffness, biased initialization
data_array_eps_truss_bias = zeros(pop_size,4,num_runs);
designs_array_eps_truss_bias = strings(pop_size,num_runs);
for i = 1:num_runs
    [data_array_eps_truss_bias(:,:,i), designs_array_eps_truss_bias(:,i)] = read_csv_file(~fibre_stiffness, 'eps_moea', feas_and_stab, biased_init, i-1);
end
[obj1_combined_eps_truss_bias, obj2_combined_eps_truss_bias, feas_combined_eps_truss_bias, stab_combined_eps_truss_bias, des_combined_eps_truss_bias] = create_combined_arrays(data_array_eps_truss_bias, designs_array_eps_truss_bias, num_runs);

% For AOS MOEA, fibre stiffness, random initialization
data_array_aos_fib_rand = zeros(pop_size,4,num_runs);
designs_array_aos_fib_rand = strings(pop_size,num_runs);
for i = 1:num_runs
    [data_array_aos_fib_rand(:,:,i), designs_array_aos_fib_rand(:,i)] = read_csv_file(fibre_stiffness, 'aos', feas_and_stab, ~biased_init, i-1);
end
[obj1_combined_aos_fib_rand, obj2_combined_aos_fib_rand, feas_combined_aos_fib_rand, stab_combined_aos_fib_rand, des_combined_aos_fib_rand] = create_combined_arrays(data_array_aos_fib_rand, designs_array_aos_fib_rand, num_runs);

% For AOS MOEA, truss stiffness, random initialization
data_array_aos_truss_rand = zeros(pop_size,4,num_runs);
designs_array_aos_truss_rand = strings(pop_size,num_runs);
for i = 1:num_runs
    [data_array_aos_truss_rand(:,:,i), designs_array_aos_truss_rand(:,i)] = read_csv_file(~fibre_stiffness, 'aos', feas_and_stab, ~biased_init, i-1);
end
[obj1_combined_aos_truss_rand, obj2_combined_aos_truss_rand, feas_combined_aos_truss_rand, stab_combined_aos_truss_rand, des_combined_aos_truss_rand] = create_combined_arrays(data_array_aos_truss_rand, designs_array_aos_truss_rand, num_runs);

% For DNF, fibre stiffness, random initialization
data_array_dnf_fib_rand = zeros(pop_size,4,num_runs);
designs_array_dnf_fib_rand = strings(pop_size,num_runs);
for i = 1:num_runs
    [data_array_dnf_fib_rand(:,:,i), designs_array_dnf_fib_rand(:,i)] = read_csv_file(fibre_stiffness, 'sc_dnf', feas_and_stab, ~biased_init, i-1);
end
[obj1_combined_dnf_fib_rand, obj2_combined_dnf_fib_rand, feas_combined_dnf_fib_rand, stab_combined_dnf_fib_rand, des_combined_dnf_fib_rand] = create_combined_arrays(data_array_dnf_fib_rand, designs_array_dnf_fib_rand, num_runs);

% For DNF, truss stiffness, random initialization
data_array_dnf_truss_rand = zeros(pop_size,4,num_runs);
designs_array_dnf_truss_rand = strings(pop_size,num_runs);
for i = 1:num_runs
    [data_array_dnf_truss_rand(:,:,i), designs_array_dnf_truss_rand(:,i)] = read_csv_file(~fibre_stiffness, 'sc_dnf', feas_and_stab, ~biased_init, i-1);
end
[obj1_combined_dnf_truss_rand, obj2_combined_dnf_truss_rand, feas_combined_dnf_truss_rand, stab_combined_dnf_truss_rand, des_combined_dnf_truss_rand] = create_combined_arrays(data_array_dnf_truss_rand, designs_array_dnf_truss_rand, num_runs);

% For ACH, fibre stiffness, random initialization
data_array_ach_fib_rand = zeros(pop_size,4,num_runs);
designs_array_ach_fib_rand = strings(pop_size,num_runs);
for i = 1:num_runs
    [data_array_ach_fib_rand(:,:,i), designs_array_ach_fib_rand(:,i)] = read_csv_file(fibre_stiffness, 'sc_ach', feas_and_stab, ~biased_init, i-1);
end
[obj1_combined_ach_fib_rand, obj2_combined_ach_fib_rand, feas_combined_ach_fib_rand, stab_combined_ach_fib_rand, des_combined_ach_fib_rand] = create_combined_arrays(data_array_ach_fib_rand, designs_array_ach_fib_rand, num_runs);

% For ACH, truss stiffness, random initialization
data_array_ach_truss_rand = zeros(pop_size,4,num_runs);
designs_array_ach_truss_rand = strings(pop_size,num_runs);
for i = 1:num_runs
    [data_array_ach_truss_rand(:,:,i), designs_array_ach_truss_rand(:,i)] = read_csv_file(~fibre_stiffness, 'sc_ach', feas_and_stab, ~biased_init, i-1);
end
[obj1_combined_ach_truss_rand, obj2_combined_ach_truss_rand, feas_combined_ach_truss_rand, stab_combined_ach_truss_rand, des_combined_ach_truss_rand] = create_combined_arrays(data_array_ach_truss_rand, designs_array_ach_truss_rand, num_runs);

%% Parameter definitions
E = 10000; % Young's Modulus for polymeric material (example: 10000 Pa)
sel = 0.05; % Unit square side length (NOT individual truss length) (example: 5 cm)
r = 50*(10^-6); % Radius for cross-sectional area of (assumed circular) truss members (example: 50 micrometers)
A = pi*(r^2); % Cross-sectional area of truss member
NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1]; 
CA_all = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 2,3; 2,4; 2,5; 2,6; 2,7; 2,8; 2,9; 3,4; 3,5; 3,6; 3,7; 3,8; 3,9; 4,5; 4,6; 4,7; 4,8; 4,9; 5,6; 5,7; 5,8; 5,9; 6,7; 6,8; 6,9; 7,8; 7,9; 8,9];
c_ratio = 1;
sidenum = 3;

%% Plot combined penalized pareto front with callback function
figure
% Plotting fibre stiffness comparison plots
plot_case = 'pen_objectives';
runcase = 0; % combined run case
f_penalized_eps_fib_rand = [obj1_combined_eps_fib_rand, obj2_combined_eps_fib_rand];
plot_pareto_seak_postproc(f_penalized_eps_fib_rand,NC,CA_all,des_combined_eps_fib_rand,feas_combined_eps_fib_rand,stab_combined_eps_fib_rand,plot_case,2,'*b');
hold on
f_penalized_eps_fib_bias = [obj1_combined_eps_fib_bias, obj2_combined_eps_fib_bias];
plot_pareto_seak_postproc(f_penalized_eps_fib_bias,NC,CA_all,des_combined_eps_fib_bias,feas_combined_eps_fib_bias,stab_combined_eps_fib_bias,plot_case,2,'ro');
hold on
f_penalized_aos_fib_rand = [obj1_combined_aos_fib_rand, obj2_combined_aos_fib_rand];
plot_pareto_seak_postproc(f_penalized_aos_fib_rand,NC,CA_all,des_combined_aos_fib_rand,feas_combined_aos_fib_rand,stab_combined_aos_fib_rand,plot_case,2,'sg');
hold on 
f_penalized_dnf_fib_rand = [obj1_combined_dnf_fib_rand, obj2_combined_dnf_fib_rand];
plot_pareto_seak_postproc(f_penalized_dnf_fib_rand,NC,CA_all,des_combined_dnf_fib_rand,feas_combined_dnf_fib_rand,stab_combined_dnf_fib_rand,plot_case,2','dm');
hold on
f_penalized_ach_fib_rand = [obj1_combined_ach_fib_rand, obj2_combined_ach_fib_rand];
plot_pareto_seak_postproc(f_penalized_ach_fib_rand,NC,CA_all,des_combined_ach_fib_rand,feas_combined_ach_fib_rand,stab_combined_ach_fib_rand,plot_case,2,'pc');
hold off
legend('eps fib rand','eps fib bias','aos fib','dnf fib','ach fib','Location','Best')

figure
% Plotting truss stiffness comparison plots
plot_case = 'pen_objectives';
runcase = 0; % combined run case
f_penalized_eps_truss_rand = [obj1_combined_eps_truss_rand, obj2_combined_eps_truss_rand];
plot_pareto_seak_postproc(f_penalized_eps_truss_rand,NC,CA_all,des_combined_eps_truss_rand,feas_combined_eps_truss_rand,stab_combined_eps_truss_rand,plot_case,2,'*b');
hold on
f_penalized_eps_truss_bias = [obj1_combined_eps_truss_bias, obj2_combined_eps_truss_bias];
plot_pareto_seak_postproc(f_penalized_eps_truss_bias,NC,CA_all,des_combined_eps_truss_bias,feas_combined_eps_truss_bias,stab_combined_eps_truss_bias,plot_case,2,'ro');
hold on
f_penalized_aos_truss_rand = [obj1_combined_aos_truss_rand, obj2_combined_aos_truss_rand];
plot_pareto_seak_postproc(f_penalized_aos_truss_rand,NC,CA_all,des_combined_aos_truss_rand,feas_combined_aos_truss_rand,stab_combined_aos_truss_rand,plot_case,2,'sg');
hold on 
f_penalized_dnf_truss_rand = [obj1_combined_dnf_truss_rand, obj2_combined_dnf_truss_rand];
plot_pareto_seak_postproc(f_penalized_dnf_truss_rand,NC,CA_all,des_combined_dnf_truss_rand,feas_combined_dnf_truss_rand,stab_combined_dnf_truss_rand,plot_case,2','dm');
hold on
f_penalized_ach_truss_rand = [obj1_combined_ach_truss_rand, obj2_combined_ach_truss_rand];
plot_pareto_seak_postproc(f_penalized_ach_truss_rand,NC,CA_all,des_combined_ach_truss_rand,feas_combined_ach_truss_rand,stab_combined_ach_truss_rand,plot_case,2,'pc');
hold off
legend('eps truss rand','eps truss bias','aos truss','dnf truss','ach truss','Location','Best')

%% Functions
function [csv_data_array, designs_array] = read_csv_file(fib_stiff, algo, feas_stab, bias_init, run_num) 
    %%%% read appropriate files 
    filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\';
    if strcmp(algo,'eps_moea')
        fileloc = 'Epsilon MOEA Runs\\';
        if fib_stiff
            filename_model = 'Fibre Stiffness code run results\\';
            if bias_init
                filename = strcat('Biased initialization\\EpsilonMOEA_emoea',num2str(run_num),'_biasedinit_fibrestiffness.csv');
            else
                filename = strcat('Random initialization\\EpsilonMOEA_emoea',num2str(run_num),'_fibrestiffness.csv');
            end
        else
            filename_model = 'Truss code run results\\';
            if bias_init
                filename = strcat('Biased initialization\\EpsilonMOEA_emoea',num2str(run_num),'_biasedinit_trussstiffness.csv');
            else
                filename = strcat('Random initialization\\EpsilonMOEA_emoea',num2str(run_num),'_trussstiffness.csv');
            end
        end
    elseif strcmp(algo,'aos') 
        switch feas_stab
            case 0
                fileloc = 'AOS MOEA Runs\\Feas and Stab False\\';
                if fib_stiff
                    filename_model = strcat('Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_fibrestiffness.csv');
                else
                    filename_model = strcat('Truss code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_trussstiffness.csv');
                end
            case 1
                fileloc = 'AOS MOEA Runs\\Feas True\\';
                if fib_stiff
                    filename_model = strcat('Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_fibrestiffness.csv');
                else
                    filename_model = strcat('Truss code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_trussstiffness.csv');
                end
            case 2
                fileloc = 'AOS MOEA Runs\\Feas and Stab True\\';
                if fib_stiff
                    filename_model = strcat('Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_fibrestiffness.csv');
                else
                    filename_model = strcat('Truss code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_trussstiffness.csv');
                end
        end
        filename = '';
    elseif strcmp(algo,'sc_dnf')
        fileloc = 'Soft Constraint Runs\\Disjunctive Normal Form Runs\\';
        if fib_stiff
            filename_model = strcat('Fibre Stiffness code run results\\Random initialization\\EpsilonMOEA_emoea_dnf',num2str(run_num),'_fibrestiffness.csv');
        else
            filename_model = strcat('Truss code run results\\Random initialization\\EpsilonMOEA_emoea_dnf',num2str(run_num),'_trussstiffness.csv');
        end
        filename = '';
    elseif strcmp(algo,'sc_ach')
        fileloc = 'Soft Constraint Runs\\Adaptive Constraint Handling Runs\\';
        if fib_stiff
            filename_model = strcat('Fibre Stiffness code run results\\Random initialization\\EpsilonMOEA_emoea_ach',num2str(run_num),'_fibrestiffness.csv');
        else
            filename_model = strcat('Truss code run results\\Random initialization\\EpsilonMOEA_emoea_ach',num2str(run_num),'_trussstiffness.csv');
        end
        filename = '';
    end
    full_filepath = strcat(filepath,fileloc,filename_model,filename);

    data_table = readtable(full_filepath,'Format','%s%f%f%f%f','HeaderLines',1);

    %%%% store retrieved data into different variables
    %%%% csv_data includes: [Pen. Obj. 1, Pen.Obj. 2, Feasibility Score,
    %%%% Stablity Score]
    pop_size =  size(data_table,1);
    csv_data = zeros(pop_size,4);
    designs = strings(pop_size);
    csv_data = data_table(:,2:end);
    designs = data_table(:,1);
    
    csv_data_array = table2array(csv_data);
    designs_array = table2array(designs);
end

function [obj1_combined, obj2_combined, feas_combined, stab_combined, designs_combined] = create_combined_arrays(data_array, designs_array, n_runs)
    n_total = 0;
    for i = 1:n_runs
        current_data_array = data_array(:,:,i);
        n_total = n_total + size(current_data_array(:,1),1);
    end
    obj1_combined = zeros(n_total,1);
    obj2_combined = zeros(n_total,1);
    feas_combined = zeros(n_total,1);
    stab_combined = zeros(n_total,1);
    designs_combined = strings(n_total,1);
    index = 1;
    for i = 1:n_runs
        current_data_array = data_array(:,:,i);
        current_designs_array = designs_array(:,i);
        n_current = size(current_data_array(:,1),1);
        obj1_combined(index:index+n_current-1,1) = current_data_array(:,1);
        obj2_combined(index:index+n_current-1,1) = current_data_array(:,2);
        feas_combined(index:index+n_current-1,1) = current_data_array(:,3);
        stab_combined(index:index+n_current-1,1) = current_data_array(:,4);
        designs_combined(index:index+n_current-1,1) = current_designs_array(:,1);
        index = index + n_current;
    end
end

function savename = get_save_name(fib_stiff,algor,run_case,bias_init,run_num)
    filename1 = 'pareto_front_';
    if fib_stiff
        filename2 = 'fibre_stiffness_';
    else
        filename2 = 'truss_stiffness_';
    end
    if strcmp(algor,'eps_moea')
        filename3 = 'eps_moea_';
    elseif strcmp(algor,'aos')
        filename3 = 'aos_';
    elseif strcmp(algor,'sc_dnf')
        filename3 = 'dnf_';
    elseif strcmp(algor,'sc_ach')
        filename3 = 'ach_';
    end
    if (run_case == 0)
        filename4 = 'combined';
    elseif (run_case == 1)
        filename4 = strcat('run',num2str(run_num));
    end
    if bias_init
        filename5 = '_biasedinit';
    else
        filename5 = '_randominit';
    end
    format = '.png';
    savename = strcat(filename1,filename2,filename3,filename4,filename5,format);
end
