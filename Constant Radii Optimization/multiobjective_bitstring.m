function [obj] = multiobjective_bitstring(x_vec,CA_all,NC,A,E,sel,r,c_ratio)
% This function computes the objectives for the bitstring design and 
% incorporates the feasibility and stability scores into the computation  
% of the objectives 
    %x_des = [0, 1, 0, 0, x_vec(1:2), 1, x_vec(3:4), 0, x_vec(5:18), ...
    %x_vec(19:end), 0, 1, 0]; % for forcing truss elements in 1,3; 1,7 and
    % 7;9
    
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
    
    %%% FEASIBILITY CONDITIONAL OBJECTIVE COMPUTATION 
    %%% Compute feasibility and stability scores for the design
    %feasibility_score = feasibility_checker_nonbinary(NC,CA_des);
    
    %sidenum = 3;
    %stability_score = stabilityTester_2D_updated(sidenum,CA_des,NC);
    
    %if (feasibility_score == 1 && stability_score == 1) 
        %%% Compute material stiffness matrix for the design
        %C = zeros(3); % units are [Pa]
        %[C_des,~,~] = generateC(sel,r,NC,CA_des,A,E,C,CA_infeasible,inf_count);
        
        % Using Fibre Stiffness Model
        %nucfac = 1;
        %[C11,C22] = fiberStiffnessModel(sel,r,E,CA_des,sidenum,nucfac);  
        %C_des = zeros(2);
        %C_des(1,1) = C11;
        %C_des(2,2) = C22;
        
        %%% Compute volume fraction for the design 
        %vf_des = calcVF(NC,CA_des,r,sel);
    
        %%% Find true objective values
        %obj = [abs((C_des(2,2)/C_des(1,1)) - c_ratio), -(C_des(2,2)/vf_des)];
        
    %elseif (feasibility_score < 1 || stability_score < 1)
        
        %%% Set objective values to penalize infeasible and unstable
        %%% designs
        %obj = -1000*[(feasibility_score),(stability_score)];
       
    %end
    
    %%% AUGMENTED LAGRANGIAN OBJECTIVE
    %alpha = 1e-4;
    %if (rem(iter,1000) == 0 && iter~= 0)
        %lambda = lambda*(1 - (alpha*iter));
    %end
    
    %disp(lambda)
    %disp(iter)  
    
    %lambda2 = 6000*(1 - (alpha*iter));
    feasibility_score = feasibility_checker_nonbinary(NC,CA_des);
    
    sidenum = 3;
    stability_score = stabilityTester_2D_updated(sidenum,CA_des,NC);
    
    %%% Compute material stiffness matrix for the design
    C = zeros(3); % units are [Pa]
    [C_des,~,~] = generateC(sel,r,NC,CA_des,A,E,C);
    
    % Using Fibre Stiffness Model
    %nucfac = 1;
    %[C11,C22] = fiberStiffnessModel(sel,r,E,CA_des,sidenum,nucfac);  
    %C_des = zeros(2);
    %C_des(1,1) = C11;
    %C_des(2,2) = C22;
    
    %%% Compute volume fraction for the design 
    vf_des = calcVF(NC,CA_des,r,sel);
    
    if (isnan(C_des(2,2)))
        C_des(2,2) = 1e-3;
        C_des(1,1) = 1e-6;
    end
    
    % objectives normalized by 15 and 6000 respectivey
    obj = [abs((C_des(2,2)/C_des(1,1)) - c_ratio)/15, -(C_des(2,2)/vf_des)/6000] - (log(abs(feasibility_score)) + log(abs(stability_score)))/2*ones(1,2);
    %obj = [abs((C_des(2,2)/C_des(1,1)) - c_ratio)/15, -(C_des(2,2)/vf_des)/6000] - (tanh(feasibility_score + 0.1) + tanh(stability_score + 0.1))/2*ones(1,2);
    %obj = [abs((C_des(2,2)/C_des(1,1)) - c_ratio)/15, -(C_des(2,2)/vf_des)/6000] - (feasibility_score + stability_score)/2*ones(1,2);
    %iter = iter + 1;
    
    %disp(iter)
    %%% ADDITIONAL FEASIBILITY+STABILITY OBJECTIVE
    %feasibility_score = feasibility_checker_nonbinary(NC,CA_des);
    
    %sidenum = 3;
    %stability_score = stabilityTester_2D_updated(sidenum,CA_des,NC);
    
    %%% Compute material stiffness matrix for the design
    %C = zeros(3); % units are [Pa]
    %[C_des,~,~] = generateC(sel,r,NC,CA_des,A,E,C,CA_infeasible,inf_count);
    
    % Using Fibre Stiffness Model
    %nucfac = 1;
    %[C11,C22] = fiberStiffnessModel(sel,r,E,CA_des,sidenum,nucfac);  
    %C_des = zeros(2);
    %C_des(1,1) = C11;
    %C_des(2,2) = C22;
    
    %%% Compute volume fraction for the design 
    %vf_des = calcVF(NC,CA_des,r,sel);
    
    %if (isnan(C_des(2,2)))
        %C_des(2,2) = 1;
        %C_des(1,1) = 1e10;
    %end
    
    %obj = [abs((C_des(2,2)/C_des(1,1)) - c_ratio), -(C_des(2,2)/vf_des), -(feasibility_score + stability_score)];
    
end

