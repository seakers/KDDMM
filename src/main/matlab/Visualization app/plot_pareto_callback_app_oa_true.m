function [] = plot_pareto_callback_app_oa_true(gcbo,eventdata,NC,CA_all,design_array,feas_array,stab_array,orient_array,f_vals_true,sel,r,E,pen_fac,sidenum,nucfac,app,filename_intpen,prob_num)
% This function is the callback for the pareto plot in the
% TrussGAVizapp app. For the selected point on the pareto front,
% the truss design plot, feasibility and stability scores, stiffness
% matrix and the true objectives are computed and the 
% Inputs: Nodal Connectivity Array NC,
%         Complete Connectivity Array CA_all,
%         Array of pareto designs design_array,
%         Array of pareto feasibility scores feas_array,
%         Array of pareto stability scores stab_array,
%         Array of pareto orientation scores orient_array
%         True objectives in the pareto front f_vals_true,
%         Number of nodes on one side sel,
%         Truss member radius r,
%         Material Young's Modulus E,
%         GUI apphandles object apphandles

axes_handle = get(gcbo,'Parent');
point_pos = get(axes_handle, 'Currentpoint');

%if prob_num=='1'
    %point_pos_comp = [point_pos(1,1),-point_pos(1,2)];
%elseif prob_num=='2'
    %point_pos_comp = [-point_pos(1,1),point_pos(1,2)];
%end

euclid_dist = zeros(size(f_vals_true,1),1);
for i = 1:size(f_vals_true,1)
    current_point = f_vals_true(i,:);
    euclid_dist(i) = sqrt((current_point(1) - point_pos(1,1))^2 + (current_point(2) - point_pos(1,2))^2);
end
[~,min_index] = min(euclid_dist);

f_true_design = f_vals_true(min_index,:);
point_design = design_array(min_index);
feas_design = feas_array(min_index);
stab_design = stab_array(min_index);
orient_design = orient_array(min_index);
point_des_char = point_design{1,1};

cla(app.ParetoAxes);
%if (prob_num == '1')
    %%f_vals_corrected = [f_vals_pen(:,1), -f_vals_pen(:,2)];
    %plot(app.ParetoAxes,f_vals_corrected(:,1),f_vals_corrected(:,2),'*','ButtonDownFcn',{@plot_pareto_callback_app_oa,NC,CA_all,design_array,feas_array,stab_array,orient_array,f_vals_pen,f_vals_true,sel,r,E,pen_fac,sidenum,nucfac,app,filename_intpen,prob_num})
    %hold(app.ParetoAxes,'on')
    %plot(app.ParetoAxes,f_pen_design(1),-f_pen_design(2),'r*')
    %hold(app.ParetoAxes,'off')
%elseif (prob_num == '2')
    %f_vals_corrected = [-f_vals_pen(:,1), f_vals_pen(:,2)];
    %plot(app.ParetoAxes,f_vals_corrected(:,1),f_vals_corrected(:,2),'*','ButtonDownFcn',{@plot_pareto_callback_app_oa,NC,CA_all,design_array,feas_array,stab_array,orient_array,f_vals_pen,f_vals_true,sel,r,E,pen_fac,sidenum,nucfac,app,filename_intpen,prob_num})
    %hold(app.ParetoAxes,'on')
    %plot(app.ParetoAxes,-f_pen_design(1),f_pen_design(2),'r*')
    %hold(app.ParetoAxes,'off')
%end
plot(app.ParetoAxes,f_vals_true(:,1),f_vals_true(:,2),'*','ButtonDownFcn',{@plot_pareto_callback_app_oa_true,NC,CA_all,design_array,feas_array,stab_array,orient_array,f_vals_true,sel,r,E,pen_fac,sidenum,nucfac,app,filename_intpen,prob_num})
hold(app.ParetoAxes,'on')
plot(app.ParetoAxes,f_true_design(1),f_true_design(2),'r*')
hold(app.ParetoAxes,'off')

design_bool = zeros(size(point_des_char,2),1);
for i = 1:size(point_des_char,2)
    design_bool(i) = str2num(point_des_char(i));
end

CA_des = CA_all(design_bool~=0,:);

visualize_truss_3x3_app(NC,CA_des,app)

%setappdata(gui_obj,'feas_score',feas_design);
%setappdata(gui_obj,'stab_score',stab_design);
%set(apphandles.hEditFeasScore,'string',num2str(feas_design));
%set(apphandles.hEditStabScore,'string',num2str(stab_design));
app.FeasibilityScoreTextArea.Value = num2str(round(feas_design,2));
app.StabilityScoreTextArea.Value = num2str(round(stab_design,2)); 
app.OrientationScoreTextArea.Value = num2str(round(orient_design,2));

A = pi*(r^2);
C = [];
[C_design,~,~] = generateC(sel,r,NC,CA_des,A,E,C);
C_design_rounded = zeros(3);
for i = 1:3
    for j = 1:3
        C_design_rounded(i,j) = round(C_design(i,j),2);
    end
end

%setappdata(gui_obj,'C_mat',C_design);
%set(apphandles.hEditCMat,'Data',C_design);
app.TrussCMat.Data = C_design_rounded;

f_true_rounded = [round(f_true_design(1),3), round(f_true_design(2),3)];

%set(apphandles.true_objs,'string',mat2str(f_true));
%setappdata(gui_obj,'true_objs',f_true);
%set(apphandles.hEditTrueObjs,'string',mat2str(f_true));
app.TrueObjectives.Data = f_true_rounded;

vf_des = calcVF(NC,CA_des,r,sel);
app.VolumeFractionTextArea.Value = num2str(round(vf_des,4));

[C11_fib, C22_fib] = fiberStiffnessModel(sel,r,E,CA_des,sidenum,nucfac);
CMat_fib = [round(C11_fib,3), round(C22_fib,3)];
app.FibreCTab.Data = CMat_fib;

end