function [state,options,optchanged] = outputfcn_savegen(options,state,flag)
% This is a custom output function that records the entire population,
% viable population, objectives and constraint scores in each generation

persistent viable_history history gen_array %current_spread previous_spread
optchanged = false;
E = 10000; % Young's Modulus for polymeric material (example: 10000 Pa)
sel = 0.05; % Unit square side length (NOT individual truss length) (example: 5 cm)
radius = 50*(10^-6); % Radius for cross-sectional area of (assumed circular) truss members (example: 50 micrometers)
A = pi*(radius^2); % Cross-sectional area of truss member
% Nodal Coordinate Vector (Standard 3x3, 2D Grid) below (each row represents a node, first column is x-coordinates, second column is y-coordinates):
NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1]; 
CA_all = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 2,3; 2,4; 2,5; 2,6; 2,7; 2,8; 2,9; 3,4; 3,5; 3,6; 3,7; 3,8; 3,9; 4,5; 4,6; 4,7; 4,8; 4,9; 5,6; 5,7; 5,8; 5,9; 6,7; 6,8; 6,9; 7,8; 7,9; 8,9];
sidenum = 3;
target_c_ratio = 1; % Ratio of C22/C11 to target
switch flag
    case 'init'
        history(:,:,1) = state.Population;
        assignin('base','gapopulationhistory',history);
        %gen_map = containers.Map('KeyType','int32','ValueType','char');
        %gen_struct = struct([]); 
        gen_array = zeros(options.PopulationSize, 4, options.MaxGenerations);
        x_curr = false(1,32);
        value_array = zeros(options.PopulationSize,4);
        count = 0;
        viable_population = false(options.PopulationSize,32); 
        for i = 1:size(state.Population,1)
            x_curr = state.Population(i,:);
            
            x_des = [x_curr(1), x_curr(2), x_curr(3), x_curr(4), x_curr(5), x_curr(6), ...
            x_curr(7), x_curr(8), x_curr(9), x_curr(10), x_curr(11), x_curr(12), ...
            x_curr(13), x_curr(14), x_curr(15), x_curr(16), x_curr(17), x_curr(3), ...
            x_curr(18), x_curr(19), x_curr(20), x_curr(21), x_curr(22), x_curr(23), ...
            x_curr(24), x_curr(25), x_curr(26), x_curr(27), x_curr(28), x_curr(29), ...
            x_curr(30), x_curr(31), x_curr(23), x_curr(1), x_curr(32), x_curr(9)]; 
            % repeatable designs
            CA_des = CA_all(x_des~=0,:);
    
            % Computing feasibility and stability scores
            feas_score = feasibility_checker_nonbinary(NC,CA_des);
            
            if (feas_score == 1)
                viable_population(count+1,:) = x_curr;
                stab_score = stabilityTester_2D_updated(sidenum, CA_des, NC);
                f_des = multiobjective_bitstring(x_curr,CA_all,NC,A,E,sel,radius,target_c_ratio);
                f_des_true = [15*(f_des(1) + log(abs(stab_score))/2), 6000*(f_des(2) + log(abs(stab_score))/2)];
                value_array(count+1,:) = string([f_des_true, feas_score, stab_score]); 
                count = count + 1;
            end
        end
        %gen_map(0) = value_array(1:count,:); 
        %gen_struct(1).gen0 = value_array(1:count,:); 
        %assignin('base','gagenstruct',gen_struct);
        gen_array(1:count,:,1) = value_array(1:count,:);
        assignin('base','gagenarray',gen_array);
        viable_history(:,:,1) = viable_population;
        assignin('base','gaviablepopulationhistory',viable_history);
        %previous_spread = state.Spread;
        
    case 'iter'
        % Update the history.
        ss = size(history,3);
        history(:,:,ss+1) = state.Population;
        assignin('base','gapopulationhistory',history);
        
        % Save the viable designs and scores from the current population
        x_curr = false(1,32);
        value_array = zeros(size(state.Population,1),4);
        count = 0;
        viable_population = strings(size(state.Population,1),32); 
        for i = 1:size(state.Population,1)
            x_curr = state.Population(i,:);
            
            x_des = [x_curr(1), x_curr(2), x_curr(3), x_curr(4), x_curr(5), x_curr(6), ...
            x_curr(7), x_curr(8), x_curr(9), x_curr(10), x_curr(11), x_curr(12), ...
            x_curr(13), x_curr(14), x_curr(15), x_curr(16), x_curr(17), x_curr(3), ...
            x_curr(18), x_curr(19), x_curr(20), x_curr(21), x_curr(22), x_curr(23), ...
            x_curr(24), x_curr(25), x_curr(26), x_curr(27), x_curr(28), x_curr(29), ...
            x_curr(30), x_curr(31), x_curr(23), x_curr(1), x_curr(32), x_curr(9)]; 
            % repeatable designs
            CA_des = CA_all(x_des~=0,:);
    
            % Computing feasibility and stability scores
            feas_score = feasibility_checker_nonbinary(NC,CA_des);
            
            if (feas_score == 1)
                viable_population(count+1,:) = x_curr;
                stab_score = stabilityTester_2D_updated(sidenum, CA_des, NC);
                f_des = multiobjective_bitstring(x_curr,CA_all,NC,A,E,sel,radius,target_c_ratio);
                f_des_true = [15*(f_des(1) + log(abs(stab_score))/2), 6000*(f_des(2) + log(abs(stab_score))/2)];
                value_array(count+1,:) = string([f_des_true, feas_score, stab_score]); 
                count = count + 1;
            end
        end
        %gen_map(state.Generation) = value_array(1:count,:); 
        %fieldname = strcat('gen',string(state.Generation));
        %gen_struct(1).fieldname = value_array(1:count,:); 
        %assignin('base','gagenstruct',gen_struct);
        gen_array(1:count,:,ss+1) = value_array(1:count,:);
        assignin('base','gagenarray',gen_array);
        viable_history(:,:,ss+1) = viable_population;
        assignin('base','gaviablepopulationhistory',viable_history);
        
        % Update spread values for termination condition
        %current_spread = state.Spread; 
        
        % Stop further computation if the change in spread is less than function
        % tolerance or maximum number of generations is reached
        %current_population = state.Population;
        %pareto_set = current_population(state.Rank == 1);
        %if ((abs(current_spread(end) - previous_spread(end)) <= options.FunctionTolerance) || (state.Generation == options.MaxGenerations))
            %state.StopFlag = 'y';
        %end
        %previous_spread = state.Spread;
        
    case 'done'
        % Include the final population in the history.
        ss = size(history,3);
        history(:,:,ss+1) = state.Population;
        assignin('base','gapopulationhistory',history);
        
        % Include final viable population in the viable history
        x_curr = false(1,32);
        value_array = strings(size(state.Population,1),4);
        count = 0;
        viable_population = false(size(state.Population,1),32); 
        for i = 1:size(state.Population,1)
            x_curr = state.Population(i,:);
            
            x_des = [x_curr(1), x_curr(2), x_curr(3), x_curr(4), x_curr(5), x_curr(6), ...
            x_curr(7), x_curr(8), x_curr(9), x_curr(10), x_curr(11), x_curr(12), ...
            x_curr(13), x_curr(14), x_curr(15), x_curr(16), x_curr(17), x_curr(3), ...
            x_curr(18), x_curr(19), x_curr(20), x_curr(21), x_curr(22), x_curr(23), ...
            x_curr(24), x_curr(25), x_curr(26), x_curr(27), x_curr(28), x_curr(29), ...
            x_curr(30), x_curr(31), x_curr(23), x_curr(1), x_curr(32), x_curr(9)]; 
            % repeatable designs
            CA_des = CA_all(x_des~=0,:);
    
            % Computing feasibility and stability scores
            feas_score = feasibility_checker_nonbinary(NC,CA_des);
            
            if (feas_score == 1)
                viable_population(count+1,:) = x_curr;
                stab_score = stabilityTester_2D_updated(sidenum, CA_des, NC);
                f_des = multiobjective_bitstring(x_curr,CA_all,NC,A,E,sel,radius,target_c_ratio);
                f_des_true = [15*(f_des(1) + log(abs(stab_score))/2), 6000*(f_des(2) + log(abs(stab_score))/2)];
                value_array(count+1,:) = string([f_des_true, feas_score, stab_score]); 
                count = count + 1;
            end
        end
        %gen_map(state.Generation) = value_array(1:count,:); 
        %fieldname = strcat('gen',string(state.Generation));
        %gen_struct(1).fieldname = value_array(1:count,:); 
        %assignin('base','gagenstruct',gen_struct);
        gen_array(1:count,:,ss+1) = value_array(1:count,:);
        gen_array_final = gen_array(:,:,state.Generation);
        assignin('base','gagenarray',gen_array_final);
        viable_history(:,:,ss+1) = viable_population;
        assignin('base','gaviablepopulationhistory',viable_history);
        
end
end