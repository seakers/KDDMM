%% Find case run with max penalized objectives
clear all
close all
clc

%% CSV Read parameters
fibre_model_used = false;
sidenum = 3;

num_runs = 30; % change based on run 
pop_size = 100; % change based on run

% partcoll_bools = [int_pen, AOS, biased_init, ACH] boolean array
% nodalprop_bools = [int_pen, AOS, biased_init, ACH] boolean array
% orient_bools = [int_pen, AOS, biased_init, ACH] boolean array

%% Search runs of a case
% Case booleans
constrad_read_case = false;
case_partcoll_bools = [false, true, true, false];
case_nodalprop_bools = [false, false, false, false];
case_orient_bools = [false, false, false, false];

n_members_total = nchoosek(sidenum^2,2); 

for i = 1:num_runs
    [data_array, ~] = read_csv_data_fullpop(fibre_model_used, case_partcoll_bools, case_nodalprop_bools, case_orient_bools, constrad_read_case, n_members_total, i-1);
    penobjs1_array = data_array(:,2);
    if ((30 < max(penobjs1_array)) && (max(penobjs1_array) < 500))
        disp(strcat('Max penalized objectives is found in run ',num2str(i-1),' of the current case'))
    end 
end

%% Functions
function [data_array, design_array] = read_csv_data_fullpop(fib_stiff, partcoll_bools, nodalprop_bools, orient_bools, constrad_read, n_total_members, run_num)
    filepath = "C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\OA run data - optimization Problem 2\\";
    methods = ["Int Pen","AOS","Bias Init","ACH"];
    
    if constrad_read
        filepath_constrad = "Constant Radii\\";
    else
        filepath_constrad = "Variable Radii\\";
    end
    
    if (partcoll_bools(2) || nodalprop_bools(2) || orient_bools(2))
        filename = "AOSMOEA_emoea_";
    else
        filename = "EpsilonMOEA_emoea_";
    end
    filepath2 = '';
    filename2 = '';
    constr_count = 0;
    for i = 1:4
        if (partcoll_bools(i) && nodalprop_bools(i) && ~orient_bools(i))
            constrained = strcat(methods(i)," - Partcoll and NodalProp\\");
            filename_snippet = strcat("pncon",num2str(i-1),"_");
        elseif (partcoll_bools(i) && ~nodalprop_bools(i) && ~orient_bools(i))
            constrained = strcat(methods(i)," - Partcoll\\");
            filename_snippet = strcat("pcon",num2str(i-1),"_");
        elseif (~partcoll_bools(i) && nodalprop_bools(i) && ~orient_bools(i))
            constrained = strcat(methods(i)," - NodalProp\\");
            filename_snippet = strcat("ncon",num2str(i-1),"_");
        elseif (~partcoll_bools(i) && ~nodalprop_bools(i) && orient_bools(i))
            constrained = strcat(methods(i)," - Orientation\\");
            filename_snippet = strcat("ocon",num2str(i-1),"_");
        elseif (partcoll_bools(i) && ~nodalprop_bools(i) && orient_bools(i))
            constrained = strcat(methods(i)," - Partcoll and Orientation\\");
            filename_snippet = strcat("pocon",num2str(i-1),"_");
        elseif (~partcoll_bools(i) && nodalprop_bools(i) && orient_bools(i))
            constrained = strcat(methods(i)," - NodalProp and Orientation\\");
            filename_snippet = strcat("nocon",num2str(i-1),"_");
        elseif (partcoll_bools(i) && nodalprop_bools(i) && orient_bools(i))
            constrained = strcat(methods(i)," - Partcoll Nodalprop and Orientation\\");
            filename_snippet = strcat("pnocon",num2str(i-1),"_");
        else
            constrained = '';
            filename_snippet = '';
            constr_count = constr_count + 1;
        end
        filepath2 = strcat(filepath2, constrained);
        filename2 = strcat(filename2, filename_snippet);
    end
    
    filepath_moea = '';
    if (constr_count == 4)
        filepath_moea = "Epsilon MOEA\\";
    end
    
    if fib_stiff
        filepath3 = "Fibre Stiffness\\";
        if constrad_read
            filename2 = strcat(filename2,"_prob2_fibre_fullpop.csv");
        else 
            filename2 = strcat(filename2,"_fibre_varrad_fullpop.csv");
        end
    else
        filepath3 = "Truss Stiffness\\";
        if constrad_read
            filename2 = strcat(filename2,"_prob2_truss_fullpop.csv");
        else
            filename2 = strcat(filename2,"_truss_varrad_fullpop.csv");
        end
    end
    
    %%%% read appropriate file 
    full_filepath = strcat(filepath,filepath_constrad,filepath2,filepath_moea,filepath3,filename,num2str(run_num),filename2);
    
    if constrad_read
        format_string = '%s';
        for j = 1:11
            format_string = strcat(format_string,'%f');
        end
        data_table = readtable(full_filepath,'Format',format_string,'HeaderLines',1);
    else
        format_string = '';
        for j = 1:(n_total_members+11)
            format_string = strcat(format_string,'%f');
        end
        data_table = readtable(full_filepath,'Format',format_string,'HeaderLines',1);
    end
    
    %%%% store retrieved data into different variables
    %%%% csv_data includes: [Pen. Obj. 1, Pen.Obj. 2, True Obj. 1, True Obj. 2, Feasibility Score,
    %%%% Connectivity Score, Stiffness Ratio Constraint, Partial Collapsibility Score, 
    %%%% Nodal Properties Score, Orientation Score]
    pop_size =  size(data_table,1);
    csv_data = zeros(pop_size,11);
    if constrad_read
        designs = strings(pop_size);
        designs = data_table(:,1);
    else
        designs = zeros(pop_size,n_total_members);
        designs = data_table(:,1:n_total_members);
    end
    csv_data = data_table(:,end-10:end);
    
    data_array = table2array(csv_data);
    design_array = table2array(designs);
end