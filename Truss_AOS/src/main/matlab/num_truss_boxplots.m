%% Aggregated Number of Trusses Comparison
clear
close all
clc

%% Read csv files
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

num_runs = 30; % change based on run 
pop_size = 100; % change based on run

% Compute volume fraction values for Epsilon MOEA runs using
% fibre stiffness model
[numtruss_struct_all_eps_fib] = get_num_trusses_vals(fibre_stiffness, epsilon_moea, feas_and_stab, biased_init, num_runs, pop_size);

% Compute volume fraction values for Epsilon MOEA runs using
% truss stiffness model
[numtruss_struct_all_eps_truss] = get_num_trusses_vals(~fibre_stiffness, epsilon_moea, feas_and_stab, biased_init, num_runs, pop_size);

% Compute volume fraction values for AOS MOEA runs using
% fibre stiffness model
[numtruss_struct_all_aos_fib] = get_num_trusses_vals(fibre_stiffness, ~epsilon_moea, feas_and_stab, biased_init, num_runs, pop_size);

% Compute volume fraction values for AOS MOEA runs using
% truss stiffness model
[numtruss_struct_all_aos_truss] = get_num_trusses_vals(~fibre_stiffness, ~epsilon_moea, feas_and_stab, biased_init, num_runs, pop_size);

%% Create boxplots for individual cases
labels = cell(num_runs,1);
for i = 1:num_runs
    labels{i,1} = num2str(i-1);
end

% volume fraction boxplots for Epsilon MOEA runs using fibre stiffness model 
[numtruss_array_all_eps_fib, numtruss_array_mean_eps_fib, numtruss_eps_fib_groups] = create_boxplot_arrays(numtruss_struct_all_eps_fib, num_runs);
figure
boxplot(numtruss_array_all_eps_fib,numtruss_eps_fib_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Number of Trusses')
title('Number of Trusses Boxplot for Epsilon MOEA runs using Fibre Stiffness Model')

% volume fraction boxplots for Epsilon MOEA runs using truss stiffness model 
[numtruss_array_all_eps_truss, numtruss_array_mean_eps_truss, numtruss_eps_truss_groups] = create_boxplot_arrays(numtruss_struct_all_eps_truss, num_runs);
figure
boxplot(numtruss_array_all_eps_truss, numtruss_eps_truss_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Number of Trusses')
title('Number of Trusses Boxplot for Epsilon MOEA runs using Truss Stiffness Model')

% volume fraction boxplots for AOS MOEA runs using fibre stiffness model 
[numtruss_array_all_aos_fib, numtruss_array_mean_aos_fib, numtruss_aos_fib_groups] = create_boxplot_arrays(numtruss_struct_all_aos_fib, num_runs);
figure
boxplot(numtruss_array_all_aos_fib, numtruss_aos_fib_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Number of Trusses')
title('Number of Trusses Boxplot for AOS MOEA runs using Fibre Stiffness Model')

% volume fraction boxplots for AOS MOEA runs using truss stiffness model 
[numtruss_array_all_aos_truss, numtruss_array_mean_aos_truss, numtruss_aos_truss_groups] = create_boxplot_arrays(numtruss_struct_all_aos_truss, num_runs);
figure
boxplot(numtruss_array_all_aos_truss, numtruss_aos_truss_groups,'Labels',labels)
xlabel('Run Number')
ylabel('Number of Trusses')
title('Number of Trusses Boxplot for AOS MOEA runs using Truss Stiffness Model')

%% Create boxpltos camparing different cases
case_labels = {'eps_fib','eps_truss','aos_fib','aos_truss'};
mean_bp_groups = [zeros(1,num_runs),ones(1,num_runs),2.*ones(1,num_runs),3.*ones(1,num_runs)];

% Plotting boxplots
numtruss_array = [numtruss_array_mean_eps_fib',numtruss_array_mean_eps_truss',numtruss_array_mean_aos_fib',numtruss_array_mean_aos_truss'];
figure 
boxplot(numtruss_array,mean_bp_groups,'Labels',case_labels);
ylabel('Number of Trusses')
title('Number of Trusses Comparison Boxplot')

%% Functions
function [num_trusses_struct_all] = get_num_trusses_vals(fib_stiff, eps_moea, const_depend, bias_init, n_runs, n_pop)
    num_trusses_struct_all = struct;
    
    %%%% read appropriate files 
    for i = 1:n_runs
        run_num = i-1;
        filepath = 'C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\';
        if eps_moea
            fileloc = 'Epsilon MOEA Runs\\';
            if fib_stiff
                filename_model = 'Fibre Stiffness code run results\\';
                if bias_init
                    filename = strcat('Biased Initialization\\EpsilonMOEA_emoea',num2str(run_num),'_biasedinit_fibrestiffness.csv');
                else
                    filename = strcat('Random Initialization\\EpsilonMOEA_emoea',num2str(run_num),'_fibrestiffness.csv');
                end
            else
                filename_model = 'Truss code run results\\';
                if bias_init
                    filename = strcat('Biased Initialization\\EpsilonMOEA_emoea',num2str(run_num),'_biasedinit_trussstiffness.csv');
                else
                    filename = strcat('Random Initialization\\EpsilonMOEA_emoea',num2str(run_num),'_trussstiffness.csv');
                end 
            end
        else 
            switch const_depend
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
    
        f_penalized = csv_data_array(:,1:2);
        pareto_bool = paretofront(f_penalized);
        designs_pareto = designs_array(pareto_bool==1);
        
        current_field = strcat('run_',num2str(i));
        n_trusses_pareto = zeros(size(designs_pareto,1),1);
        
        for i = 1:size(designs_pareto,1)
            current_design = designs_pareto{i};
            design_bool = zeros(size(current_design,2),1);
            for j = 1:size(design_bool,1)
                design_bool(j) = str2num(current_design(j));
            end
            n_trusses_pareto(i,1) = nnz(design_bool);
        end
        
        num_trusses_struct_all.(current_field) = n_trusses_pareto;
        
    end
end

function [val_array, mean_val_array, bp_groups] = create_boxplot_arrays(val_structs, n_runs)
    val_array = [];
    bp_groups = [];
    mean_val_array = zeros(n_runs,1);
    group_count = 0;
    for i = 1:n_runs
        current_field = strcat('run_',num2str(i));
        current_val_array = val_structs.(current_field);
        val_array = [val_array;current_val_array];
        bp_groups = [bp_groups;group_count.*ones(size(current_val_array,1),1)];
        mean_val_array(i) = mean(current_val_array);
        group_count = group_count + 1;
    end
end









