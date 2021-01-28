% CONNECTIVITY CONSTRAINT--NPBC
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
function conConstScore = connectivityConstraint_NPBC_2D(sidenum,NC,CA,sel,biasFac)
    % Initialize variables
    conConstScore = 1; 
    ND = NC./sel;
    
    % Add up counters based on nodal connectivities
    [N,~] = histcounts(CA,size(NC,1));
    
    % Loop through each node
    for i = 1:1:size(NC,1)
        % Isolate elements originating/ending at a given node
        indione = CA(:,1) == i; inditwo = CA(:,2) == i;
        mCAone = CA(indione,:); mCAtwo = CA(inditwo,:);
        mCA = ...
            [setdiff(mCAone,mCAtwo,'rows');setdiff(mCAtwo,mCAone,'rows')];
        
        % Consider elements present due to repeatability
        if (ND(i,1) == 0) && (ND(i,2) == 0) % Node in bottom left corner
            % Identify opposite node
            oppnode = i+((sidenum^2)-1);
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
                
            % Identify adjacent nodes
            badjnode = i+(sidenum-1);
            aadjnode = i+(sidenum*(sidenum-1));
            
            % Find elements connecting to/from adjacent nodes
            bindone = CA(:,1) == badjnode; bindtwo = CA(:,2) == badjnode;
            bnCAone = CA(bindone,:); bnCAtwo = CA(bindtwo,:);
            bnCA = [setdiff(bnCAone,bnCAtwo,'rows');...
                    setdiff(bnCAtwo,bnCAone,'rows')];
            aindone = CA(:,1) == aadjnode; aindtwo = CA(:,2) == aadjnode;
            anCAone = CA(aindone,:); anCAtwo = CA(aindtwo,:);
            anCA = [setdiff(anCAone,anCAtwo,'rows');...
                    setdiff(anCAtwo,anCAone,'rows')];
                
            % Based on location of node i relative to adjacent nodes, 
            % eliminate duplicate elements
            for q = 1:1:size(bnCA,1)
                if (ND(bnCA(q,1),1) ~= ND(bnCA(q,2),1)) && ...
                   (ND(bnCA(q,1),2) ~= ND(bnCA(q,2),2))
                    mCA = [mCA;bnCA(q,:)];
                    N(i) = N(i) + 1;
                end
            end
            for q = 1:1:size(anCA,1)
                if (ND(anCA(q,1),1) ~= ND(anCA(q,2),1)) && ...
                   (ND(anCA(q,1),2) ~= ND(anCA(q,2),2))
                    mCA = [mCA;anCA(q,:)];
                    N(i) = N(i) + 1;
                end
            end
            
            % Add connections from opposite node to mCA
            mCA = [mCA;onCA];      
            N(i) = N(i) + size(onCA,1); 
            
        elseif (ND(i,1) == 0) && (ND(i,2) == 1) % Node in top left corner
            % Identify opposite node
            oppnode = i+((sidenum-1)^2); 
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
                
            % Identify adjacent nodes
            badjnode = i-(sidenum-1);
            aadjnode = i+(sidenum*(sidenum-1));
            
            % Find elements connecting to/from adjacent nodes
            bindone = CA(:,1) == badjnode; bindtwo = CA(:,2) == badjnode;
            bnCAone = CA(bindone,:); bnCAtwo = CA(bindtwo,:);
            bnCA = [setdiff(bnCAone,bnCAtwo,'rows');...
                    setdiff(bnCAtwo,bnCAone,'rows')];
            aindone = CA(:,1) == aadjnode; aindtwo = CA(:,2) == aadjnode;
            anCAone = CA(aindone,:); anCAtwo = CA(aindtwo,:);
            anCA = [setdiff(anCAone,anCAtwo,'rows');...
                    setdiff(anCAtwo,anCAone,'rows')];
                
            % Based on location of node i relative to adjacent nodes, 
            % eliminate duplicate elements
            for q = 1:1:size(bnCA,1)
                if (ND(bnCA(q,1),1) ~= ND(bnCA(q,2),1)) && ...
                   (ND(bnCA(q,1),2) ~= ND(bnCA(q,2),2))
                    mCA = [mCA;bnCA(q,:)];
                    N(i) = N(i) + 1;
                end
            end
            for q = 1:1:size(anCA,1)
                if (ND(anCA(q,1),1) ~= ND(anCA(q,2),1)) && ...
                   (ND(anCA(q,1),2) ~= ND(anCA(q,2),2))
                    mCA = [mCA;anCA(q,:)];
                    N(i) = N(i) + 1;
                end
            end
            
            % Add connections from opposite node to mCA
            mCA = [mCA;onCA];      
            N(i) = N(i) + size(onCA,1); 
            
        elseif (ND(i,1) == 1) && (ND(i,2) == 0) % Node in bot. right corner
            % Identify opposite node
            oppnode = i-((sidenum-1)^2); 
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
            
            % Identify adjacent nodes
            badjnode = i+(sidenum-1);
            aadjnode = i-(sidenum*(sidenum-1));
            
            % Find elements connecting to/from adjacent nodes
            bindone = CA(:,1) == badjnode; bindtwo = CA(:,2) == badjnode;
            bnCAone = CA(bindone,:); bnCAtwo = CA(bindtwo,:);
            bnCA = [setdiff(bnCAone,bnCAtwo,'rows');...
                    setdiff(bnCAtwo,bnCAone,'rows')];
            aindone = CA(:,1) == aadjnode; aindtwo = CA(:,2) == aadjnode;
            anCAone = CA(aindone,:); anCAtwo = CA(aindtwo,:);
            anCA = [setdiff(anCAone,anCAtwo,'rows');...
                    setdiff(anCAtwo,anCAone,'rows')];
                
            % Based on location of node i relative to adjacent nodes, 
            % eliminate duplicate elements
            for q = 1:1:size(bnCA,1)
                if (ND(bnCA(q,1),1) ~= ND(bnCA(q,2),1)) && ...
                   (ND(bnCA(q,1),2) ~= ND(bnCA(q,2),2))
                    mCA = [mCA;bnCA(q,:)];
                    N(i) = N(i) + 1;
                end
            end
            for q = 1:1:size(anCA,1)
                if (ND(anCA(q,1),1) ~= ND(anCA(q,2),1)) && ...
                   (ND(anCA(q,1),2) ~= ND(anCA(q,2),2))
                    mCA = [mCA;anCA(q,:)];
                    N(i) = N(i) + 1;
                end
            end
            
            % Add connections from opposite node to mCA
            mCA = [mCA;onCA];      
            N(i) = N(i) + size(onCA,1); 
            
        elseif (ND(i,1) == 1) && (ND(i,2) == 1) % Node in top right corner
            % Identify opposite node
            oppnode = i-((sidenum^2)-1);
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
            
            % Identify adjacent nodes
            badjnode = i-(sidenum-1);
            aadjnode = i-(sidenum*(sidenum-1));
            
            % Find elements connecting to/from adjacent nodes
            bindone = CA(:,1) == badjnode; bindtwo = CA(:,2) == badjnode;
            bnCAone = CA(bindone,:); bnCAtwo = CA(bindtwo,:);
            bnCA = [setdiff(bnCAone,bnCAtwo,'rows');...
                    setdiff(bnCAtwo,bnCAone,'rows')];
            aindone = CA(:,1) == aadjnode; aindtwo = CA(:,2) == aadjnode;
            anCAone = CA(aindone,:); anCAtwo = CA(aindtwo,:);
            anCA = [setdiff(anCAone,anCAtwo,'rows');...
                    setdiff(anCAtwo,anCAone,'rows')];
                
            % Based on location of node i relative to adjacent nodes, 
            % eliminate duplicate elements
            for q = 1:1:size(bnCA,1)
                if (ND(bnCA(q,1),1) ~= ND(bnCA(q,2),1)) && ...
                   (ND(bnCA(q,1),2) ~= ND(bnCA(q,2),2))
                    mCA = [mCA;bnCA(q,:)];
                    N(i) = N(i) + 1;
                end
            end
            for q = 1:1:size(anCA,1)
                if (ND(anCA(q,1),1) ~= ND(anCA(q,2),1)) && ...
                   (ND(anCA(q,1),2) ~= ND(anCA(q,2),2))
                    mCA = [mCA;anCA(q,:)];
                    N(i) = N(i) + 1;
                end
            end
            
            % Add connections from opposite node to mCA
            mCA = [mCA;onCA];      
            N(i) = N(i) + size(onCA,1);  
            
        elseif (ND(i,1) == 0) % Edge node is on left side
            % Identify opposite node
            oppnode = i+(sidenum*(sidenum-1));
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
        
            % Based on location of node i/oppnode, eliminate duplicate
            % elements
            for q = 1:1:size(onCA,1)
                if ND(onCA(q,1),1) ~= ND(onCA(q,2),1)
                    mCA = [mCA;onCA(q,:)];
                    N(i) = N(i) + 1;
                end
            end
        elseif (ND(i,1) == 1) % Edge node is on right side
            % Identify opposite node
            oppnode = i-(sidenum*(sidenum-1));
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
        
            % Based on location of node i/oppnode, eliminate duplicate
            % elements
            for q = 1:1:size(onCA,1)
                if ND(onCA(q,1),1) ~= ND(onCA(q,2),1)
                    mCA = [mCA;onCA(q,:)];
                    N(i) = N(i) + 1;
                end
            end
        elseif (ND(i,2) == 0) % Edge node is on bottom side
            % Identify opposite node
            oppnode = i+(sidenum-1);
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
        
            % Based on location of node i/oppnode, eliminate duplicate
            % elements
            for q = 1:1:size(onCA,1)
                if ND(onCA(q,1),2) ~= ND(onCA(q,2),2)
                    mCA = [mCA;onCA(q,:)];
                    N(i) = N(i) + 1;
                end
            end
        elseif (ND(i,2) == 1) % Edge node is on top side
            % Identify opposite node
            oppnode = i-(sidenum-1);
            
            % Find elements connecting to/from oppnode
            indone = CA(:,1) == oppnode; indtwo = CA(:,2) == oppnode;
            onCAone = CA(indone,:); onCAtwo = CA(indtwo,:);
            onCA = [setdiff(onCAone,onCAtwo,'rows');...
                    setdiff(onCAtwo,onCAone,'rows')];
        
            % Based on location of node i/oppnode, eliminate duplicate
            % elements
            for q = 1:1:size(onCA,1)
                if ND(onCA(q,1),2) ~= ND(onCA(q,2),2)
                    mCA = [mCA;onCA(q,:)];
                    N(i) = N(i) + 1;
                end
            end
        end
        
        % Find sums of element components in x,y of elements originating
        % from the current node
        xsum = 0; ysum = 0;
        for j = 1:1:size(mCA,1)
            x1 = NC(mCA(j,1),1); x2 = NC(mCA(j,2),1);
            y1 = NC(mCA(j,1),2); y2 = NC(mCA(j,2),2);
            xsum = xsum + abs(x2-x1);
            ysum = ysum + abs(y2-y1);
        end
        
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