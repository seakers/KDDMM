function [] = plot_pareto_callback_gui(gcbo,eventdata,NC,CA_all,design_array,feas_array,stab_array,f_vals,sel,r,E,pen_fac,handles,gui_obj)
% This function is the callback for the pareto plot in the
% TrussGAVisualizer_simple gui. For the selected point on the pareto front,
% the truss design plot, feasibility and stability scores, stiffness
% matrix and the true objectives are computed and the 
% Inputs: Nodal Connectivity Array NC,
%         Complete Connectivity Array CA_all,
%         Array of pareto designs design_array,
%         Array of pareto feasibility scores feas_array,
%         Array of pareto stability scores stab_array,
%         Penalized objectives in the pareto front f_vals,
%         Number of nodes on one side sel,
%         Truss member radius r,
%         Material Young's Modulus E,
%         GUI handles object handles

point_pos = get(gca, 'Currentpoint');

point_pos_comp = [point_pos(1,1),-point_pos(1,2)];

euclid_dist = zeros(size(f_vals,1),1);
for i = 1:size(f_vals,1)
    current_point = f_vals(i,:);
    euclid_dist(i) = sqrt((current_point(1) - point_pos_comp(1))^2 + (current_point(2) - point_pos_comp(2))^2);
end
[~,min_index] = min(euclid_dist);

f_pen_design = f_vals(min_index,:);
point_design = design_array(min_index);
feas_design = feas_array(min_index);
stab_design = stab_array(min_index);
point_des_char = point_design{1,1};

cla(handles.hTrussFig);
f_vals_true = [f_vals(:,1), -f_vals(:,2)];
plot(f_vals_true(:,1),f_vals_true(:,2),'*','Parent',handles.hParetoFig,'ButtonDownFcn',{@plot_pareto_callback_gui,NC,CA_all,design_array,feas_array,stab_array,f_vals,sel,r,E,pen_fac,handles,gui_obj})
hold on
plot(f_pen_design(1),-f_pen_design(2),'r*','Parent',handles.hParetoFig)
hold off

design_bool = zeros(size(point_des_char,2),1);
for i = 1:size(point_des_char,2)
    design_bool(i) = str2num(point_des_char(i));
end

CA_des = CA_all(design_bool~=0,:);

visualize_truss_3x3_gui(NC,CA_des,handles.hTrussFig)

%set(handles.feas_score,'string',num2str(feas_design));
%set(handles.stab_score,'string',num2str(stab_design));
setappdata(gui_obj,'feas_score',feas_design);
setappdata(gui_obj,'stab_score',stab_design);
set(handles.hEditFeasScore,'string',num2str(feas_design));
set(handles.hEditStabScore,'string',num2str(stab_design));

A = pi*(r^2);
C = [];
[C_design,~,~] = generateC(sel,r,NC,CA_des,A,E,C);

%set(handles.C_mat,'string',mat2str(C_design));
setappdata(gui_obj,'C_mat',C_design);
set(handles.hEditCMat,'Data',C_design);

penalty = (log10(abs(feas_design)) + log10(abs(stab_design)))/2;
f_true = [15*(f_pen_design(1) + pen_fac*penalty), -8500*(f_pen_design(2) + pen_fac*penalty)];

%set(handles.true_objs,'string',mat2str(f_true));
setappdata(gui_obj,'true_objs',f_true);
set(handles.hEditTrueObjs,'string',mat2str(f_true));

end