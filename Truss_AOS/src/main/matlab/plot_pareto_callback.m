function [] = plot_pareto_callback(gcbo, eventdata, NC, CA_all, design_array, feas_array, stab_array, f_vals)
% This is the graph callback function for the pareto front plot

point_pos = get(gca, 'Currentpoint');

point_pos_comp = [point_pos(1,1),-point_pos(1,2)];


%f_obj1 = f_vals(:,1);
%f_obj2 = f_vals(:,2);

%k1 = find(f_obj1==point_pos_comp(1));
%k2 = find(f_obj2==point_pos_comp(2));

%k_similar = intersect(k1,k2);

%point_design = design_array(k_similar);
%feas_design = feas_array(k_similar);
%stab_design = stab_array(k_similar);

euclid_dist = zeros(size(f_vals,1),1);
for i = 1:size(f_vals,1)
    current_point = f_vals(i,:);
    euclid_dist(i) = sqrt((current_point(1) - point_pos_comp(1))^2 + (current_point(2) - point_pos_comp(2))^2);
end
[~,min_index] = min(euclid_dist);

point_design = design_array(min_index);
feas_design = feas_array(min_index);
stab_design = stab_array(min_index);
point_des_char = point_design{1,1};

design_bool = zeros(size(point_des_char,2),1);
for i = 1:size(point_des_char,2)
    design_bool(i) = str2num(point_des_char(i));
end

CA_des = CA_all(design_bool~=0,:);

visualize_truss_3x3(NC,CA_des)

disp(strcat('Feasibility Score = ',num2str(feas_design)))
disp(strcat('Stability Score = ',num2str(stab_design)))
disp(strcat('Connectivity Array = ',mat2str(CA_des)))
fprintf('\n')

end

