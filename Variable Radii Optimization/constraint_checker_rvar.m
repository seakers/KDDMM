function [c,c_eq] = constraint_checker_rvar(x_rvar,nucFac,sel,CA_all)
% This function is the nonlinear constraint function for the real coded
% variable radii GA optimization.

    % Calculated Inputs
    sidenum = (2*nucFac) + 1; % Number of nodes on each size of square grid
    
    % Compute full design vector
    x_rvar_full = [x_rvar(1:8), x_rvar(2), x_rvar(9:16), x_rvar(11), x_rvar(1), x_rvar(4)];
    
    x_binary = (x_rvar_full >= 1e-6); % binary thresholded design vector for
    % feasibility and stability score computation
    
    CA_des = CA_all(x_binary~=0,:);

    % Generate vector with nodal coordinates
    NC = generateNC(sel,sidenum);
    
    % Compute feasibility and stability scores
    feasibility_score = feasibility_checker_nonbinary(NC,CA_des);
    stability_score = stabilityTester_2D_updated(sidenum,CA_des,NC);
    
    % Setting equality constraint values
    c = [];
    %c_eq = (log(abs(feasibility_score)) + log(abs(stability_score)))/2*ones(1,2); 
    c_eq = [log(abs(feasibility_score)), log(abs(stability_score))];

end