function [feasibility_score, feasibility_bool] = feasibility_checker_boolean(x,NC,CA_all)
% This function is the boolean version of the feasiblity_checker_nonbinary
% function
% Inputs: design x
%         nodal position matrix NC
%         full connectivity array CA_all
    % Converting design vector to CA matrix
    
    % x is a double vector, so converting to logical vector
    %x_vec = x>0.5;
    %k = find(x_vec);
    
    %x_vec = [0, 1, 0, 0, x(1:2), 1, x(3:4), 0, x(5:18), x(19:end), 0, 1, 0]; % for forcing truss elements in 1,3; 1,7 and
    %7;9
    %x_vec = x; % general case
    %x_vec = [x(1), 0, x(2:3), 0, 0, 0, 0, x(4:7), 0, 0, 0, 0, x(8:9), 0, 0, 0, ...
        %x(10), 0, x(11:12), 0, x(13:16), 0, x(17:19), 0, x(20)]; % only adjacent node connections allowed
    
    %CA_des = CA_all(x_vec~=0,:);
    
    x_des = [x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10), ...
        x(11), x(12), x(13), x(14), x(15), x(16), x(17), x(3), x(18), x(19), ...
        x(20), x(21), x(22), x(23), x(24), x(25), x(26), x(27), x(28), x(29), ...
        x(30), x(31), x(23), x(1), x(32), x(9)]; 
    % repeatable designs  
    
    CA_des = CA_all(x_des~=0,:);
    
    %disp(CA_des)
    %n_truss_des = size(CA_des,1);
    %disp(n_truss_des)
    
    % Feasibility check
    feasibility_score = feasibility_checker_nonbinary(NC,CA_des);
    
    % Assign value to constraint boolean
    feasibility_bool = true;
    if feasibility_score < 1
        feasibility_bool = false;
    end
    
end

