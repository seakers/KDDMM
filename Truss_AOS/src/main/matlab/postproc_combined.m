%% Plot combined pareto front
clear all
close all
clc

%% Read results from csv files
%%%% set to: 
%       true - to read fibre stiffness evaluation run results
%       false - to read truss stiffness evaluation run results
fibre_stiffness = false;  

%%%% set to: 
%       true - to read Epsilon MOEA run results
%       false - to read AOS MOEA run results
epsilon_moea = false;

%%%% set to: 
%       0 - both feasibility and stability is not used in both operators
%       1 - only feasibility is used ine AddTruss operators
%       2 - both feasibility and stability are used in both operators
feas_and_stab = 0;

num_runs = 30; % change based on run 
pop_size = 100; % change based on run
data_array = zeros(pop_size,4,num_runs);
designs_array = strings(pop_size,num_runs);

for i = 1:num_runs
    [data_array(:,:,i), designs_array(:,i)] = read_csv_file(fibre_stiffness, epsilon_moea, feas_and_stab, i-1);
end

%% Form combined objectives and scores vectors
n_total = 0;
for i = 1:num_runs
    current_data_array = data_array(:,:,i);
    n_total = n_total + size(current_data_array(:,1),1);
end
obj1_combined = zeros(n_total,1);
obj2_combined = zeros(n_total,1);
feas_combined = zeros(n_total,1);
stab_combined = zeros(n_total,1);
designs_combined = strings(n_total,1);
index = 1;
for i = 1:num_runs
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
f_penalized = [obj1_combined, obj2_combined];
plot_case = 'pen_objectives';
runcase = 0; % combined run case
figure
plot_pareto_seak_postproc(f_penalized,NC,CA_all,designs_combined,feas_combined,stab_combined,plot_case,2);
%saveas(gcf,get_save_name(fibre_stiffness,epsilon_moea,runcase,0));

%% Plot individual pareto fronts with callbcak function
plot_case = 'pen_objectives';
runcase = 1; % individual run case
for i = 1:num_runs
    obj1_current = data_array(:,1,i);
    obj2_current = data_array(:,2,i);
    feas_current = data_array(:,3,i);
    stab_current = data_array(:,4,i);
    designs_current = designs_array(:,i);
    f_penalized_current = [obj1_current, obj2_current];
    figure
    plot_pareto_seak_postproc(f_penalized_current,NC,CA_all,designs_current,feas_current,stab_current,plot_case,2);
    %saveas(gcf,get_save_name(fibre_stiffness,epsilon_moea,runcase,i-1));
end

%% Functions
function [csv_data_array, designs_array] = read_csv_file(fib_stiff, eps_moea, feas_stab, run_num) 
    %%%% read appropriate files 
    filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\';
    if eps_moea
        fileloc = 'Epsilon MOEA Runs\\';
        if fib_stiff
            filename = strcat('Fibre Stiffness code run results\\EpsilonMOEA_emoea',num2str(run_num),'_fibrestiffness.csv');
        else
            filename = strcat('Truss code run results\\EpsilonMOEA_emoea',num2str(run_num),'_trussstiffness.csv');
        end
    else 
        switch feas_stab
            case 0
                fileloc = 'AOS MOEA Runs\\Feas and Stab False\\';
                if fib_stiff
                    filename = strcat('Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_fibrestiffness.csv');
                    %credit_filename = 'Fibre Stiffness code run results\\constraint_adaptive0.credit';
                    %quality_filename = 'Fibre Stiffness code run results\\constraint_adaptive0.qual';
                else
                    filename = strcat('Truss code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_trussstiffness.csv');
                    %credit_filename = 'Truss code run results\\constraint_adaptive0.credit';
                    %quality_filename = 'Truss code run results\\constraint_adaptive0.qual';
                end
            case 1
                fileloc = 'AOS MOEA Runs\\Feas True\\';
                if fib_stiff
                    filename = strcat('Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_fibrestiffness.csv');
                    %credit_filename = 'Fibre Stiffness code run results\\constraint_adaptive0.credit';
                    %quality_filename = 'Fibre Stiffness code run results\\constraint_adaptive0.qual';
                else
                    filename = strcat('Truss code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_trussstiffness.csv');
                    %credit_filename = 'Truss code run results\\constraint_adaptive0.credit';
                    %quality_filename = 'Truss code run results\\constraint_adaptive0.qual';
                end
            case 2
                fileloc = 'AOS MOEA Runs\\Feas and Stab True\\';
                if fib_stiff
                    filename = strcat('Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_fibrestiffness.csv');
                    %credit_filename = 'Fibre Stiffness code run results\\constraint_adaptive0.credit';
                    %quality_filename = 'Fibre Stiffness code run results\\constraint_adaptive0.qual';
                else
                    filename = strcat('Truss code run results\\AOSMOEA_constraint_adaptive',num2str(run_num),'_trussstiffness.csv');
                    %credit_filename = 'Truss code run results\\constraint_adaptive0.credit';
                    %quality_filename = 'Truss code run results\\constraint_adaptive0.qual';
                end
        end
    end
    full_filepath = strcat(filepath,fileloc,filename);
    %full_filepath_credit = strcat(filepath,fileloc,credit_filename);
    %full_filepath_quality = strcat(filepath,fileloc,quality_filename);

    data_table = readtable(full_filepath,'Format','%s%f%f%f%f','HeaderLines',1);
    %credit_fileId = fopen(full_filepath_credit,'r');
    %quality_fileId = fopen(full_filepath_quality,'r');
    %array_size = [6 Inf];
    %credits = fscanf(credit_fileId,'%c',array_size);
    %quality = fscanf(quality_fileId,'%c',array_size);

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

function savename = get_save_name(fib_stiff,eps_moea,run_case,run_num)
    filename1 = 'pareto_front_';
    if fib_stiff
        filename2 = 'fibre_stiffness_';
    else
        filename2 = 'truss_stiffness_';
    end
    if eps_moea
        filename3 = 'eps_moea_';
    else
        filename3 = 'aos_';
    end
    if (run_case == 0)
        filename4 = 'combined';
    elseif (run_case == 1)
        filename4 = strcat('run',num2str(run_num));
    end
    format = '.png';
    savename = strcat(filename1,filename2,filename3,filename4,format);
end

