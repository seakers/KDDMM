function [stable] = stabilityTester_2D_boolean(sidenum,x_vec,CA_all,NC)
% This function assesses whether a given design is stable or not
% Inputs: sidenum - number of nodes on each side 
%         x - design vector
%         NC - nodeal connectivity matrix
% Output: stable - bollean stability
%x_des = [0, 1, 0, 0, x_vec(1:2), 1, x_vec(3:4), 0, x_vec(5:18), ...
%x_vec(19:end), 0, 1, 0]; % for forcing truss elements in 1,3; 1,7 and 7;9
    
%x_des = x_vec; % general case
    
%x_des = [x_vec(1), 0, x_vec(2:3), 0, 0, 0, 0, x_vec(4:7), 0, 0, 0, 0, ... 
%x_vec(8:9), 0, 0, 0, x_vec(10), 0, x_vec(11:12), 0, x_vec(13:16), ...
%0, x_vec(17:19), 0, x_vec(20)]; % only adjacent node connections allowed

x_des = [x_vec(1), x_vec(2), x_vec(3), x_vec(4), x_vec(5), x_vec(6), ...
        x_vec(7), x_vec(8), x_vec(9), x_vec(10), x_vec(11), x_vec(12), ...
        x_vec(13), x_vec(14), x_vec(15), x_vec(16), x_vec(17), x_vec(3), ...
        x_vec(18), x_vec(19), x_vec(20), x_vec(21), x_vec(22), x_vec(23), ...
        x_vec(24), x_vec(25), x_vec(26), x_vec(27), x_vec(28), x_vec(29), ...
        x_vec(30), x_vec(31), x_vec(23), x_vec(1), x_vec(32), x_vec(9)]; 
        % repeatable designs  

CA_des = CA_all(x_des~=0,:);

stable = false;

stability_score = stabilityTester_2D_updated(sidenum,CA_des,NC);
if stability_score == 1
    stable = true;

end

