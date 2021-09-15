% Constraint Checker - Non binary constraints
function [c,c_eq] = constraint_checker_nonbinary(x,NC,CA_all,inf_count,CA_infeasible)
% This function is intended to check that a given truss design is
% legitimate (feasible, stable, and having the minimum number of truss elements)
    % Converting design vector to CA matrix
    % x is a double vector, so converting to logical vector
    x_vec = x>0.5;
    %x_vec = x >= 5; %double vector, IntCon 
    %k = find(x_vec);
    CA_des = CA_all(x_vec~=0,:);
    %disp(CA_des)
    %n_truss_des = size(CA_des,1);
    %disp(n_truss_des)
    
    % Feasibility check
    feasibility_score = feasibility_checker_nonbinary(NC,CA_des);

    % Stability check
    sidenum = 3;
    stability_score = stabilityTester_2D_updated(sidenum,CA_des,NC);

    % Truss element number check
    %b = 10; % minimum number of truss elements
    %truss_num_score = max([0, 1/(n_truss_des - b)]);
    
    % Inequality constraint
    c = [(-feasibility_score),(-stability_score)];
    %c = [(-feasibility_score),(-stability_score), truss_num_score];
    %c = [];

    % Dummy Equality Constraint
    %c_eq = [feasibility_score-1,stability_score-1];
    c_eq = [];

    %if ~all(c)
        %disp('Infeasible CA')
        %%disp(CA_des)
        %disp('Feasibility')
        %disp(feasibility_score)
        %disp('Stability')
        %disp(stability_score)
        %%disp('Number of truss elements scores')
        %%disp(truss_num_score)
        %inf_count = inf_count + 1;
        %CA_infeasible{inf_count} = CA_des;
    %end
    
end