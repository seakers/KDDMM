function [obj] = multiobjective_rvar(x_rvar,nucFac,sel,E,CA_all,c_ratio)
% This function computes the objectives for the GA to optimize. This is for 
% the real-coded case where the design vector consists of the radius values
% for the trusses between the two corresponding nodes

    % Calculated Inputs
    sidenum = (2*nucFac) + 1; % Number of nodes on each size of square grid
    
    % Compute full design vector
    x_rvar_full = [x_rvar(1:8), x_rvar(2), x_rvar(9:16), x_rvar(11), x_rvar(1), x_rvar(4)];
    
    x_binary = (x_rvar_full >= 1e-6); % binary thresholded design vector for
    % feasibility and stability score computation
    x_rvar_updated = x_rvar_full.*x_binary; % setting the radii below the 
    % threshold to zero
    
    CA_des = CA_all(x_binary~=0,:);
    
    x_rvar_des = nonzeros(x_rvar_updated);
    Avar = pi.*(x_rvar_des.^2); % Cross-sectional areas of truss members
    
    % Generate vector with nodal coordinates
    NC = generateNC(sel,sidenum);
    
    % Compute feasibility and stability scores
    feasibility_score = feasibility_checker_nonbinary(NC,CA_des);
    stability_score = stabilityTester_2D_updated(sidenum,CA_des,NC);
    
    % Develop C-matrix from K-matrix (functions below)
    C = [];
    [C_des,~,~] = generateC(sel,x_rvar_des,NC,CA_des,Avar,E,C,sidenum);
    C_des = C_des./nucFac;

    if (isnan(C_des(2,2)))
        C_des(2,2) = 1e-3;
        C_des(1,1) = 1e-6;
    end
    
    % Calculate volume fraction
    volFrac = calcVF(NC,CA_des,x_rvar_des,sel,sidenum);
    
    % Compute objectives with interior penalties
    obj = [abs((C_des(2,2)/C_des(1,1)) - c_ratio)/15, -(C_des(2,2)/volFrac)/6000] - (log(abs(feasibility_score)) + log(abs(stability_score)))/2*ones(1,2);
    
end









