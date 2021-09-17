function [obj] = multiobjective_evMOGA(des,param)
    %%% Extracting parameters from param
    CA_all = param.CA_all;
    NC = param.NC;
    A = param.A;
    E = param.E;
    sel = param.sel;
    r = param.r;
    c_ratio = param.c_ratio;
    CA_infeasible = param.CA_infeasible;
    inf_count = param.inf_count;
    
    %%% Converting design vector to CA matrix
    x_vec = des>0.5;
    %k = find(des);
    CA_des = CA_all(x_vec~=0,:);
    
    %%% AUGMENTED LAGRANGIAN OBJECTIVE
    lambda = 100 ;
    feasibility_score = feasibility_checker_nonbinary(NC,CA_des);
    
    sidenum = 3;
    stability_score = stabilityTester_2D_updated(sidenum,CA_des,NC);
    
    %%% Compute material stiffness matrix for the design
    C = zeros(3); % units are [Pa]
    [C_des,~,~] = generateC(sel,r,NC,CA_des,A,E,C,CA_infeasible,inf_count);
        
    %%% Compute volume fraction for the design 
    vf_des = calcVF(NC,CA_des,r,sel);
    
    if (isnan(C_des(2,2)))
        C_des(2,2) = 1;
        C_des(1,1) = 1e10;
    end
    
    obj = [abs((C_des(2,2)/C_des(1,1)) - c_ratio), -(C_des(2,2)/vf_des)] - lambda*[(feasibility_score + stability_score), (feasibility_score + stability_score)];
        
end