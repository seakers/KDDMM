function [Population] = customCreateFunction(GenomeLength, FitnessFcn, options)
% This function is a custom binary population creation function for the
% multiobjective GA. It adds four feasible and stable designs and generates
% the rest randomly.
% Inputs: GenomeLength: number of independent variables for the fitness
% function
%         FitnessFcn: function evaluating the objectives for each design
%         options: Additional options

n_pop = options.PopulationSize;
Population = false(n_pop, GenomeLength);

% All except least four designs are randomly generated
for i = 1:n_pop-4
    Population(i,:) = randi([0 1], GenomeLength, 1);
end

% Last four designs are feasible and stable  
Population(end-3,:) = [1,0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1];
Population(end-2,:) = [1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,0,0,0,1,1,0,1];
Population(end-1,:) = [1,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1];
Population(end,:) = [1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,1,0,0];
    
end

