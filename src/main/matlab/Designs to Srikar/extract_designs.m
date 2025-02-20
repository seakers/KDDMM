%% Extract and determine unqie designs with full feasibility and stability from csv files
clear
clc

%% Extract fully feasible and stable designs from csv files
filename1 = "EpsilonMOEA_emoea_";
filename2 = "fscon0__fibre.csv";

n_runs = 30;

designs_feas_and_stab_allruns = [];

for i = 1:n_runs
    full_filepath = strcat(filename1,num2str(i-1),filename2);
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
    
    feas_array = csv_data_array(:,3);
    stab_array = csv_data_array(:,4);
    
    designs_array_feas_and_stab = designs_array((feas_array == 1) & (stab_array == 1));
    designs_feas_and_stab_allruns = [designs_feas_and_stab_allruns;designs_array_feas_and_stab];
end

%% Remove non-unique designs 
designs_allruns_unique = unique(designs_feas_and_stab_allruns);

%% Convert to Connectivity Array 
% Boolean string representation of design is converted to CA matrix and
% stored in a structure
CA_all = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 2,3; 2,4; 2,5; 2,6; 2,7; 2,8; 2,9; 3,4; 3,5; 3,6; 3,7; 3,8; 3,9; 4,5; 4,6; 4,7; 4,8; 4,9; 5,6; 5,7; 5,8; 5,9; 6,7; 6,8; 6,9; 7,8; 7,9; 8,9];
CA_allruns_struct = struct; 

for i = 1:size(designs_allruns_unique)
    current_des_boolstring = designs_allruns_unique{i};
    design_bool = zeros(size(current_des_boolstring,2),1);
    for j = 1:size(current_des_boolstring,2)
        design_bool(j) = str2num(current_des_boolstring(j));
    end

    CA_des = CA_all(design_bool~=0,:);
    current_field = strcat('design',num2str(i));
    CA_allruns_struct.(current_field) = CA_des;
end
    
