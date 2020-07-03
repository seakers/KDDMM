%%% Testing the feasiblity and stability of neighbouring designs of
%%% feasible and stable designs 
clear
close all
clc

%% Defining the problem constants
E = 10000; % Young's Modulus for polymeric material (example: 10000 Pa)
sel = 0.05; % Unit square side length (NOT individual truss length) (example: 5 cm)
r = 50*(10^-6); % Radius for cross-sectional area of (assumed circular) truss members (example: 50 micrometers)
A = pi*(r^2); % Cross-sectional area of truss member
% Nodal Coordinate Vector (Standard 3x3, 2D Grid) below (each row represents a node, first column is x-coordinates, second column is y-coordinates):
NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1]; 
CA_all = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 2,3; 2,4; 2,5; 2,6; 2,7; 2,8; 2,9; 3,4; 3,5; 3,6; 3,7; 3,8; 3,9; 4,5; 4,6; 4,7; 4,8; 4,9; 5,6; 5,7; 5,8; 5,9; 6,7; 6,8; 6,9; 7,8; 7,9; 8,9];
sidenum = 3;

%% The feasible and stable architectures
population_feas(1,:) = [1,0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1];
population_feas(2,:) = [1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,0,0,0,1,1,0,1];
population_feas(3,:) = [1,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1];
population_feas(4,:) = [1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,1,0,0];

%% Checking the neighbouring designs
feas_and_stab = false(4,size(population_feas(1,:),2));
n_feas_and_stab = zeros(4,1);
for i = 1:4
    des = population_feas(i,:);
    for j = 1:36
        des_mod = des;
        des_mod(j) = ~des(j);
        [~,feas_des_mod] = feasibility_checker_boolean(des_mod,NC,CA_all);
        stab_des_mod = stabilityTester_2D_boolean(sidenum,des_mod,CA_all,NC);
        if (feas_des_mod && stab_des_mod)
            feas_and_stab(i,j) = true;
        end
    end
    n_feas_and_stab(i,1) = nnz(feas_and_stab(i,:));
end

      
