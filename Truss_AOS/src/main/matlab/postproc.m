%% Epsilon MOEA post proc
clear all
close all
clc

%% Read results csv file
%%%% set to: 
%       true - to read fibre stiffness evaluation run results
%       false - to read truss stiffness evaluation run results
fibre_stiffness = true;  

%%%% set to: 
%       true - to read Epsilon MOEA run results
%       false - to read AOS MOEA run results
epsilon_moea = true;

%%%% set to: 
%       0 - both feasibility and stability is not used in both operators
%       1 - only feasibility is used ine AddTruss operators
%       2 - both feasibility and stability are used in both operators
feas_and_stab = 0;

%%%% set to:
%       true - to read biased initialization results
%       false - to read random initialization results
biased_init = true;

%%%% read appropriate files 
filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\';
if epsilon_moea
    fileloc = 'Epsilon MOEA Runs\\';
    if fibre_stiffness
        filename_model = 'Fibre Stiffness code run results\\';
        if biased_init
            filename = 'Biased initialization\\EpsilonMOEA_emoea0_biasedinit_fibrestiffness.csv';
        else
            filename = 'Random initialization\\EpsilonMOEA_emoea0_fibrestiffness.csv';
        end
    else
        filename_model = 'Truss code run results\\';
        if biased_init
            filename = 'Biased initialization\\EpsilonMOEA_emoea0_biasedinit_fibrestiffness.csv';
        else
            filename = 'Random initialization\\EpsilonMOEA_emoea0_fibrestiffness.csv';
        end
    end
else 
    switch feas_and_stab
        case 0
            fileloc = 'AOS MOEA Runs\\Feas and Stab False\\';
            if fibre_stiffness
                filename_model = 'Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive0_fibrestiffness.csv';
            else
                filename_model = 'Truss code run results\\AOSMOEA_constraint_adaptive0_trussstiffness.csv';
            end
        case 1
            fileloc = 'AOS MOEA Runs\\Feas True\\';
            if fibre_stiffness
                filename_model = 'Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive0_fibrestiffness.csv';
                
            else
                filename_model = 'Truss code run results\\AOSMOEA_constraint_adaptive0_trussstiffness.csv';
            end
        case 2
            fileloc = 'AOS MOEA Runs\\Feas and Stab True\\';
            if fibre_stiffness
                filename_model = 'Fibre Stiffness code run results\\AOSMOEA_constraint_adaptive0_fibrestiffness.csv';
            else
                filename_model = 'Truss code run results\\AOSMOEA_constraint_adaptive0_trussstiffness.csv';
            end
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

%% Compute the true objective values 
designs_array = designs{:,:};
f_penalized = csv_data{:,1:2};
feas_scores = csv_data{:,3}; 
stab_scores = csv_data{:,4};
f_true = zeros(pop_size, 2);
pen_fac = 1;
if fibre_stiffness
    pen_fac = 1.5;
end
for i = 1:pop_size
    penalty = (log10(abs(feas_scores(i))) + log10(abs(stab_scores(i))))/2;
    f_true(i,:) = [15*(f_penalized(i,1) + pen_fac*penalty), 85000*(f_penalized(i,2) + pen_fac*penalty)];
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

%% Plot pareto front using penalized objectives with callback function
plot_case = 'pen_objectives';
figure
plot_pareto_seak_postproc(f_penalized,NC,CA_all,designs_array,feas_scores,stab_scores,plot_case,2);

%% Plot Pareto Front using true objectives with callback function
% isolate only fully feasible designs
%fully_feas = feas_scores == 1;
%f_true_feas = f_true(fully_feas ~= 0,:);
%design_feas = designs_array(fully_feas ~= 0);
%feas_scores_feas = feas_scores(fully_feas ~= 0);
%stab_scores_feas = stab_scores(fully_feas ~= 0);

% use all designs 
f_true_feas = f_true;
design_feas = designs_array;
feas_scores_feas = feas_scores;
stab_scores_feas = stab_scores;
plot_case = 'true_objectives';

figure
plot_pareto_seak_postproc(f_true_feas,NC,CA_all,design_feas,feas_scores_feas,stab_scores_feas,plot_case,2);

%% Compare evaluations of fibre stiffness results with truss code
if fibre_stiffness
    obj_true = zeros(pop_size,2);
    for i = 1:pop_size
        current_des = designs_array(i);
        current_des_char = current_des{1,1}; 
        design_bool = zeros(size(current_des_char,2),1);
        for j = 1:size(current_des_char,2)
            design_bool(j) = str2num(current_des_char(j));
        end
        CA_des = CA_all(design_bool~=0,:);    
        C = zeros(3); % units are [Pa]
        [C_des,~,~] = generateC(sel,r,NC,CA_des,A,E,C);
        vf_des = calcVF(NC,CA_des,r,sel);
        if (isnan(C_des(2,2)))
            C_des(2,2) = 1e-3;
            C_des(1,1) = 1e-6;
        end
        obj_true(i,:) = [abs((C_des(2,2)/C_des(1,1)) - c_ratio), C_des(2,2)/vf_des];
    end
end

