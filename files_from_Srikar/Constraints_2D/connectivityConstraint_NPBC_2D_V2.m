% CONNECTIVITY CONSTRAINT--NPBC--Truss Model-Specific
% This function determines whether or not designs have sufficient
% connections at each node to ensure stability as interpreted by models.  
% Connectivity at edge/corner nodes accounts for unit cell repetition.  
% This function is intended for use with non-periodic boundary conditions
% -------------------------------------------------------------------------
% Designs are scored for each violation.  Because designs on smaller node
% grids are prone to having inherently higher scores, a bias factor can be
% applied to scores when used with such designs.  The scoring variable, 
% conConstScore, is on [0,1] in increments of 0.1
% -------------------------------------------------------------------------
function conConstScore = connectivityConstraint_NPBC_2D_V2(NC,CA,biasFac)
    % Initialize variables
    conConstScore = 1; 
    
    % Add up counters based on nodal connectivities
    [N,~] = histcounts(CA,size(NC,1));
    
    % Loop through each node
    for i = 1:1:size(NC,1)
        % Determine whether node has sufficient connectivity components
        if (N(i) < 2)
            % Node is an unstable connection point
            conConstScore = conConstScore - 0.1;
            if conConstScore < 0.1
                return
            end
        end 
    end
    
    % Account for bias factor
    conConstScore = conConstScore.*biasFac;
end