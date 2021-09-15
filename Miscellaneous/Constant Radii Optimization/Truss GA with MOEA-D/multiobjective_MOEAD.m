function [obj] = multiobjective_MOEAD(des,CA_all,NC,A,E,sel,r,c_ratio,CA_infeasible,inf_count)
    % Converting design vector to CA matrix
    x_vec = des>0.5;
    %x_vec = des>=5; % double vector, IntCon
    %k = find(des);
    CA_des = CA_all(x_vec~=0,:);
    
    %%% AUGMENTED LAGRANGIAN OBJECTIVE
    lambda = 100 ;
    feasibility_score = feasibility_checker_nonbinary(NC,CA_des);
    
    sidenum = 3;
    stability_score = stabilityTester_2D_updated(sidenum,CA_des,NC);
    
    C = zeros(3); % units are [Pa]
    %print('Calculating C matrix for current design......')
    [C_des,~,~] = generateC(sel,r,NC,CA_des,A,E,C,CA_infeasible,inf_count);
    
    if (isnan(C_des(2,2)))
        C_des(2,2) = 1;
        C_des(1,1) = 1e10;
    end
    
    % Using Fibre Stiffness Model
    %nucfac = 1;
    %sidenum = 3;
    %[C11,C22] = fiberStiffnessModel(sel,r,E,CA_des,sidenum,nucfac);  
    %C_des = zeros(2);
    %C_des(1,1) = C11;
    %C_des(2,2) = C22;
    
    %print('Calculating Volume Fraction for current design......')
    vf_des = calcVF(NC,CA_des,r,sel);
    obj = [abs((C_des(2,2)/C_des(1,1)) - c_ratio); -(abs(C_des(2,2)/vf_des))] - lambda*[(feasibility_score + stability_score); (feasibility_score + stability_score)];;
end