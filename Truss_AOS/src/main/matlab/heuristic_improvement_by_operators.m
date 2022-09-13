%% Check heuristic satisfaction after application of heuristic operators
clear
close all
clc

%% Extract and analyze data for requisite problem
truss_problem = false; % true -> truss problem, false -> artery problem

filepath = "C:\\SEAK Lab\\SEAK Lab Github\\KD3M3\\Truss_AOS\\result\\";

filename = "operator_heuristic_satisfaction";
if truss_problem
    filename = strcat(filename,"_truss.csv");
else
    filename = strcat(filename,"_artery.csv");
end
filepath = strcat(filepath,filename);

format = '%s%s%f%f%s%f%f%s%f%f%s%f%f'; 
% [Full_design_initial, Full_design_partcoll, partcoll_old, partcoll_new, Full_design_nodalprop, nodalprop_old, nodalprop_new, 
% Full_design_orient, orient_old, orient_new, Full_design_inters, inters_old, inters_new]

data_table = readtable(filepath,'Format',format,'HeaderLines',1);

partcoll_old = table2array(data_table(:,3));
partcoll_new = table2array(data_table(:,4));
nodalprop_old = table2array(data_table(:,6));
nodalprop_new = table2array(data_table(:,7));
orient_old = table2array(data_table(:,9));
orient_new = table2array(data_table(:,10));
inters_old = table2array(data_table(:,12));
inters_new = table2array(data_table(:,13));

% Peform Welch's t-test test to compare old and new datasets for different heuristics 
[~, p_partcoll] = ttest2(partcoll_old, partcoll_new, 'VarType', 'unequal', 'Tail', 'left');
[~, p_nodalprop] = ttest2(nodalprop_old, nodalprop_new, 'VarType', 'unequal', 'Tail', 'left');
[~, p_orient] = ttest2(orient_old, orient_new, 'VarType', 'unequal', 'Tail', 'left');
[~, p_inters] = ttest2(inters_old, inters_new, 'VarType', 'unequal', 'Tail', 'left');
    
% Plot boxplots
labels_before = repmat({'Before Operator'},size(partcoll_old,1),1);
labels_after = repmat({'After Operator'},size(partcoll_new,1),1);
figure
subplot(2,2,1)
boxplot([partcoll_old; partcoll_new],[labels_before; labels_after])
title('Partial Collapsibility')
subplot(2,2,2)
boxplot([nodalprop_old; nodalprop_new],[labels_before; labels_after])
title('Nodal Properties')
subplot(2,2,3)
boxplot([orient_old; orient_new],[labels_before; labels_after])
title('Orientation')
subplot(2,2,4)
boxplot([inters_old; inters_new],[labels_before; labels_after])
title('Intersection')
sgtitle('Heuristic Improvement by Operators','FontSize',10) 

%% Compute I1 using heuristic operators
I1_partcoll = mean(partcoll_new - partcoll_old);
I1_nodalprop = mean(nodalprop_new - nodalprop_old);
I1_orient = mean(orient_new - orient_old);
I1_inters = mean(inters_new - inters_old);

%% Compute Cohen's d for heuristics
n_partcoll_old = size(partcoll_old, 1);
n_partcoll_new = size(partcoll_new, 1);
n_nodalprop_old = size(nodalprop_old, 1);
n_nodalprop_new = size(nodalprop_new, 1);
n_orient_old = size(orient_old, 1);
n_orient_new = size(orient_new, 1);
n_inters_old = size(inters_old, 1);
n_inters_new = size(inters_new, 1);

s_pooled_partcoll = sqrt(((n_partcoll_old - 1)*var(partcoll_old) + (n_partcoll_new - 1)*var(partcoll_new))/(n_partcoll_old + n_partcoll_new - 2));
s_pooled_nodalprop = sqrt(((n_nodalprop_old - 1)*var(nodalprop_old) + (n_nodalprop_new - 1)*var(nodalprop_new))/(n_nodalprop_old + n_nodalprop_new - 2));
s_pooled_orient = sqrt(((n_orient_old - 1)*var(orient_old) + (n_orient_new - 1)*var(orient_new))/(n_orient_old + n_orient_new - 2));
s_pooled_inters = sqrt(((n_inters_old - 1)*var(inters_old) + (n_inters_new - 1)*var(inters_new))/(n_inters_old + n_inters_new - 2));

d_partcoll = (mean(partcoll_new) - mean(partcoll_old))/s_pooled_partcoll;
d_nodalprop = (mean(nodalprop_new) - mean(nodalprop_old))/s_pooled_nodalprop;
d_orient = (mean(orient_new) - mean(orient_old))/s_pooled_orient;
d_inters = (mean(inters_new) - mean(inters_old))/s_pooled_inters;