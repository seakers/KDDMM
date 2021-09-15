% PARTIAL COLLAPSIBILITY HEURISTIC
% This function determines whether or not a design is prone to partial
% collapse under shear loading
% -------------------------------------------------------------------------
% Designs are scored for each violation.  Because designs on smaller node
% grids are prone to having inherently higher scores, a bias factor can be
% applied to scores when used with such designs.  The scoring variable, 
% pcScore, is on [0,1] in increments of 0.1
% -------------------------------------------------------------------------
function pcScore = partCollapseHeuristic_2D(sidenum,CA,NC,sel,biasFac)
    % Initialize variables
    pcScore = 1; 
    ND = NC./sel;

    % Iterate through slices in x-direction
    for ix = 1:1:(sidenum-1)
        % Identify nodes on the surface of the current slice
        lowerbound = (ix-1)*(1/(sidenum-1));
        upperbound = (ix)*(1/(sidenum-1));
        nodeslower = find(ND(:,1) == lowerbound); 
        nodesupper = find(ND(:,1) == upperbound); 
        allnodes = [nodeslower;nodesupper];

        % Isolate all elements connecting to/from all these nodes
        mCAone = CA(ismember(CA(:,1),allnodes),:);
        mCAtwo = CA(ismember(CA(:,2),allnodes),:);
        mCA = intersect(mCAone,mCAtwo,'rows');

        % Also include diagonal elements starting on the slice but ending
        % elsewhere on the grid
        for ib = 1:1:size(nodeslower,1) % For nodeslower
            tCAone = CA(ismember(CA(:,1),nodeslower(ib)),:);
            for io = 1:1:size(tCAone,1)
                if ND(tCAone(io,2),1) >= ((2/(sidenum-1))+ND(ib,1))
                    mCA = [mCA;tCAone(io,:)];
                end
            end
            tCAtwo = CA(ismember(CA(:,2),nodeslower(ib)),:);
            for io = 1:1:size(tCAtwo,1)
                if ND(tCAtwo(io,1),1) >= ((2/(sidenum-1))+ND(ib,1))
                    mCA = [mCA;tCAtwo(io,:)];
                end
            end
        end
        for ib = 1:1:size(nodesupper,1) % For nodesupper
            tCAone = CA(ismember(CA(:,1),nodesupper(ib)),:);
            for io = 1:1:size(tCAone,1)
                if ND(tCAone(io,2),1) <= (ND(ib,1)-(2/(sidenum-1)))
                    mCA = [mCA;tCAone(io,:)];
                end
            end
            tCAtwo = CA(ismember(CA(:,2),nodesupper(ib)),:);
            for io = 1:1:size(tCAtwo,1)
                if ND(tCAtwo(io,1),1) <= (ND(ib,1)-(2/(sidenum-1)))
                    mCA = [mCA;tCAtwo(io,:)];
                end
            end
        end

        % Isolate and test diagonal elements for feasibility
        xyblocker = false;
        for j = 1:1:size(mCA,1)
            x1 = NC(mCA(j,1),1); x2 = NC(mCA(j,2),1);
            y1 = NC(mCA(j,1),2); y2 = NC(mCA(j,2),2);
            if (abs(x2-x1) > 0) && (abs(y2-y1) > 0)
                xyblocker = true;
            end
        end

        % Check logical operator for stability pass/fail
        if xyblocker == false
            pcScore = pcScore - (biasFac*0.1);
            if pcScore < 0.1
                return
            end
        end
    end

    % Iterate through slices in y-direction
    for iy = 1:1:(sidenum-1)
        % Identify nodes on the surface of the current slice
        lowerbound = (iy-1)*(1/(sidenum-1));
        upperbound = (iy)*(1/(sidenum-1));
        nodeslower = find(ND(:,2) == lowerbound); 
        nodesupper = find(ND(:,2) == upperbound); 
        allnodes = [nodeslower;nodesupper];

        % Isolate all elements connecting to/from all these nodes
        mCAone = CA(ismember(CA(:,1),allnodes),:);
        mCAtwo = CA(ismember(CA(:,2),allnodes),:);
        mCA = intersect(mCAone,mCAtwo,'rows');

        % Also include diagonal elements starting on the slice but ending
        % elsewhere on the grid
        for ib = 1:1:size(nodeslower,1) % For nodeslower
            tCAone = CA(ismember(CA(:,1),nodeslower(ib)),:);
            for io = 1:1:size(tCAone,1)
                if ND(tCAone(io,2),2) >= ((2/(sidenum-1))+ND(ib,2))
                    mCA = [mCA;tCAone(io,:)];
                end
            end
            tCAtwo = CA(ismember(CA(:,2),nodeslower(ib)),:);
            for io = 1:1:size(tCAtwo,1)
                if ND(tCAtwo(io,1),2) >= ((2/(sidenum-1))+ND(ib,2))
                    mCA = [mCA;tCAtwo(io,:)];
                end
            end
        end
        for ib = 1:1:size(nodesupper,1) % For nodesupper
            tCAone = CA(ismember(CA(:,1),nodesupper(ib)),:);
            for io = 1:1:size(tCAone,1)
                if ND(tCAone(io,2),2) <= (ND(ib,2)-(2/(sidenum-1)))
                    mCA = [mCA;tCAone(io,:)];
                end
            end
            tCAtwo = CA(ismember(CA(:,2),nodesupper(ib)),:);
            for io = 1:1:size(tCAtwo,1)
                if ND(tCAtwo(io,1),2) <= (ND(ib,2)-(2/(sidenum-1)))
                    mCA = [mCA;tCAtwo(io,:)];
                end
            end
        end

        % Isolate and test diagonal elements for feasibility
        yxblocker = false;
        for j = 1:1:size(mCA,1)
            x1 = NC(mCA(j,1),1); x2 = NC(mCA(j,2),1);
            y1 = NC(mCA(j,1),2); y2 = NC(mCA(j,2),2);
            if (abs(x2-x1) > 0) && (abs(y2-y1) > 0)
                yxblocker = true;
            end
        end

        % Check logical operator for stability pass/fail
        if yxblocker == false
            pcScore = pcScore - (biasFac*0.1);
            if pcScore < 0.1
                return
            end
        end
    end
end

