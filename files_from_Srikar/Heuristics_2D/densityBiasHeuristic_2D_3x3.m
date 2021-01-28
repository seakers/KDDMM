% Function to calculate density bias score
% -----------------------------------------------------------------
% Sample values to test code (copy-paste into command window):
%{
clc;    
close all; 
clear;
scaleFac = 1;
sel = 0.05; 
% Case 1, 1 unit cell (3x3 grid)
CA = [1,2;2,3;1,4;1,5;2,5;3,5;3,6;4,5;5,6;4,7;5,7;5,8;5,9;6,9;7,8;8,9];
rvar = (50*(10^-6)).*ones(1,size(CA,1));
%}
%{
ISSUES:
(1) vF values for all 4 sections are unequal, likely due to different area
vectors
%}
function dbScore = densityBiasHeuristic_2D_3x3(CA,rvar,sel,scaleFac)
    % Inputs
    Avar = pi.*(rvar.^2); % Cross-sectional areas of truss members
    sidenum = 3; smalldim = 2;
    
    % Generate vector with nodal coordinates
    NC = generateNC(sel,sidenum);
    
    % Modify Avar for edge members for all members
    Avar = modifyAreas(Avar,CA,NC,sidenum);
    
    % Sort CA by first column values, and Avar & rvar based off of CA
    [CA,idx] = sortrows(CA); Avar = Avar(idx); rvar = rvar(idx);

    % Identify members in left and right half of design
    [lCA,lAvar,lrvar,rCA,rAvar,rrvar] = ...
        splitDesignLR(NC,CA,sidenum,rvar,Avar,sel);
    
    % Find volume fraction for each half
    lvF = calcVF(NC,lCA,lrvar,lAvar,sel,sidenum,smalldim,1,1);
    rvF = calcVF(NC,rCA,rrvar,rAvar,sel,sidenum,smalldim,1,2);
    
    % Identify members in top and bottom half of design
    [tCA,tAvar,trvar,bCA,bAvar,brvar] = ...
            splitDesignTB(NC,CA,sidenum,rvar,Avar,sel);
        
    % Find volume fraction for each half
    tvF = calcVF(NC,tCA,trvar',tAvar',sel,sidenum,smalldim,2,2);
    bvF = calcVF(NC,bCA,brvar',bAvar',sel,sidenum,smalldim,2,1);
    
    % Find differences between left & right and top & bottom VFs, and then
    % find average difference
    lrdiff = abs(lvF-rvF); tbdiff = abs(tvF-bvF); 
    avgdiff = (lrdiff+tbdiff)./2;
    
    % Score the average difference on e^-x
    dbScore = (exp(-(scaleFac*avgdiff)));
    Avars = {lAvar,rAvar,bAvar,tAvar};
    rvars = {lrvar,rrvar,brvar,trvar};
end

%-------%
% FUNCTION TO GENERATE NODAL COORDINATES BASED ON GRID SIZE
function NC = generateNC(sel,sidenum)
    notchvec = linspace(0,1,sidenum);
    NC = [];
    for i = 1:1:sidenum
        for j = 1:1:sidenum
            NC = [NC;notchvec(i),notchvec(j)];
        end
    end
    NC = sel.*NC;
end

% FUNCTION TO MODIFY AREAS FOR EDGE MEMBERS
function Avar = modifyAreas(Avar,CA,NC,sidenum)
    % Identify edge nodes
    edgenodes = [(1:1:sidenum),((2*sidenum):sidenum:...
                 ((sidenum^2)-sidenum)),((sidenum+1):sidenum:...
                 ((sidenum^2)-(2*sidenum)+1)),(((sidenum^2)-sidenum+1):...
                 1:(sidenum^2))];
             
    % Identify members connecting solely to edge nodes
    edgeconn1 = ismember(CA(:,1),edgenodes);
    edgeconn2 = ismember(CA(:,2),edgenodes);
    edgeconnectors = (edgeconn1 & edgeconn2);
    
    % Isolate edge members based on angle
    edgelogical = [edgeconnectors,edgeconnectors];
    CAedgenodes = CA.*edgelogical;
    CAedgenodes = CAedgenodes(any(CAedgenodes,2),:);
    x1 = NC(CAedgenodes(:,1),1); x2 = NC(CAedgenodes(:,2),1);
    y1 = NC(CAedgenodes(:,1),2); y2 = NC(CAedgenodes(:,2),2);
    L = sqrt(((x2-x1).^2)+((y2-y1).^2));
    angles = rad2deg(abs(acos((x2-x1)./L)));
    CAedgy = [];
    for i = 1:1:size(CAedgenodes,1)
        if (angles(i) == 0) || (angles(i) == 90)
            CAedgy = [CAedgy;CAedgenodes(i,:)];
        end
    end
    
    % Find and modify areas belonging to edge members
    if isempty(CAedgy)
        % Do nothing
    else
        edgemembers = ismember(CA,CAedgy,'rows');
        selectAreas = Avar'.*edgemembers;
        k = find(selectAreas);
        Avar(k) = Avar(k)./2;
    end
end

% FUNCTION TO SPLIT DESIGN IN HALF, VERTICALLY
function [lCA,lAvar,lrvar,rCA,rAvar,rrvar] = ...
            splitDesignLR(NC,CA,sidenum,rvar,Avar,sel)
    % Initialize variables
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
        
        % Isolate other properties of elements in this slice
        if ix == 1
            lCA = mCA;
            lAvar = nonzeros(Avar.*ismember(CA,mCA,'rows')');
            lrvar = nonzeros(rvar.*ismember(CA,mCA,'rows')');
            
            % Halve areas of members at inner edge of slice
            nuCAone = CA(ismember(CA(:,1),nodesupper),:);
            nuCAtwo = CA(ismember(CA(:,2),nodesupper),:);
            nuCA = intersect(nuCAone,nuCAtwo,'rows');
            chopindices = find(ismember(mCA,nuCA,'rows'));
            lAvar(chopindices) = lAvar(chopindices)./2;
            
        elseif ix == 2
            rCA = mCA;
            rAvar = nonzeros(Avar.*ismember(CA,mCA,'rows')');
            rrvar = nonzeros(rvar.*ismember(CA,mCA,'rows')');
            
            % Halve areas of members at inner edge of slice
            nlCAone = CA(ismember(CA(:,1),nodeslower),:);
            nlCAtwo = CA(ismember(CA(:,2),nodeslower),:);
            nlCA = intersect(nlCAone,nlCAtwo,'rows');
            chopindices = find(ismember(mCA,nlCA,'rows'));
            rAvar(chopindices) = rAvar(chopindices)./2;
            
        end
        % Also include diagonal elements starting on the slice but ending
        % elsewhere on the grid (assigning them half the area)
        for ib = 1:1:size(nodeslower,1) % For nodeslower
            tCAone = CA(ismember(CA(:,1),nodeslower(ib)),:);
            for io = 1:1:size(tCAone,1)
                if ND(tCAone(io,2),1) > (0.5*ix)
                    if ix == 1
                        lCA = [lCA;tCAone(io,:)];
                        lAvar = [lAvar;...
                            (Avar(ismember(CA,tCAone(io,:),'rows'))./2)];
                        lrvar = [lrvar;...
                            (rvar(ismember(CA,tCAone(io,:),'rows'))./2)];
                    elseif ix == 2
                        rCA = [rCA;tCAone(io,:)];
                        rAvar = [rAvar;...
                            (Avar(ismember(CA,tCAone(io,:),'rows'))./2)];
                        rrvar = [rrvar;...
                            (rvar(ismember(CA,tCAone(io,:),'rows'))./2)];
                    end
                end
            end
            tCAtwo = CA(ismember(CA(:,2),nodeslower(ib)),:);
            for io = 1:1:size(tCAtwo,1)
                if ND(tCAtwo(io,1),1) > (0.5*ix)
                    if ix == 1
                        lCA = [lCA;tCAtwo(io,:)];
                        lAvar = [lAvar;...
                            (Avar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                        lrvar = [lrvar;...
                            (rvar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                    elseif ix == 2
                        rCA = [rCA;tCAtwo(io,:)];
                        rAvar = [rAvar;...
                            (Avar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                        rrvar = [rrvar;...
                            (rvar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                    end
                end
            end
        end
        for ib = 1:1:size(nodesupper,1) % For nodesupper
            tCAone = CA(ismember(CA(:,1),nodesupper(ib)),:);
            for io = 1:1:size(tCAone,1)
                if ND(tCAone(io,2),1) < (0.5*ix)
                    if ix == 1
                        lCA = [lCA;tCAone(io,:)];
                        lAvar = [lAvar;...
                            (Avar(ismember(CA,tCAone(io,:),'rows'))./2)];
                        lrvar = [lrvar;...
                            (rvar(ismember(CA,tCAone(io,:),'rows'))./2)];
                    elseif ix == 2
                        rCA = [rCA;tCAone(io,:)];
                        rAvar = [rAvar;...
                            (Avar(ismember(CA,tCAone(io,:),'rows'))./2)];
                        rrvar = [rrvar;...
                            (rvar(ismember(CA,tCAone(io,:),'rows'))./2)];
                    end
                end
            end
            tCAtwo = CA(ismember(CA(:,2),nodesupper(ib)),:);
            for io = 1:1:size(tCAtwo,1)
                if ND(tCAtwo(io,1),1) < (0.5*ix)
                    if ix == 1
                        lCA = [lCA;tCAtwo(io,:)];
                        lAvar = [lAvar;...
                            (Avar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                        lrvar = [lrvar;...
                            (rvar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                    elseif ix == 2
                        rCA = [rCA;tCAtwo(io,:)];
                        rAvar = [rAvar;...
                            (Avar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                        rrvar = [rrvar;...
                            (rvar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                    end
                end
            end
        end
    end
    [rCA,ria,~] = unique(rCA,'rows');
    rAvar = rAvar(ria); rrvar = rrvar(ria);
    [lCA,lia,~] = unique(lCA,'rows');
    lAvar = lAvar(lia); lrvar = lrvar(lia);
end

% FUNCTION TO SPLIT DESIGN IN HALF, HORIZONTALLY
function [tCA,tAvar,trvar,bCA,bAvar,brvar] = ...
            splitDesignTB(NC,CA,sidenum,rvar,Avar,sel)
    % Initialize variables
    ND = NC./sel;
        
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
        
        % Isolate other properties of elements in this slice
        if iy == 1
            bCA = mCA;
            bAvar = nonzeros(Avar.*ismember(CA,mCA,'rows')');
            brvar = nonzeros(rvar.*ismember(CA,mCA,'rows')');
            
            % Halve areas of members at inner edge of slice
            nuCAone = CA(ismember(CA(:,1),nodesupper),:);
            nuCAtwo = CA(ismember(CA(:,2),nodesupper),:);
            nuCA = intersect(nuCAone,nuCAtwo,'rows');
            chopindices = find(ismember(mCA,nuCA,'rows'));
            bAvar(chopindices) = bAvar(chopindices)./2;
            
        elseif iy == 2
            tCA = mCA;
            tAvar = nonzeros(Avar.*ismember(CA,mCA,'rows')');
            trvar = nonzeros(rvar.*ismember(CA,mCA,'rows')');
            
            % Halve areas of members at inner edge of slice
            nlCAone = CA(ismember(CA(:,1),nodeslower),:);
            nlCAtwo = CA(ismember(CA(:,2),nodeslower),:);
            nlCA = intersect(nlCAone,nlCAtwo,'rows');
            chopindices = find(ismember(mCA,nlCA,'rows'));
            tAvar(chopindices) = tAvar(chopindices)./2;
            
        end
        %{
        % Also include diagonal elements starting on the slice but ending
        % elsewhere on the grid (assigning them half the area)
        for ib = 1:1:size(nodeslower,1) % For nodeslower
            tCAone = CA(ismember(CA(:,1),nodeslower(ib)),:);
            for io = 1:1:size(tCAone,1)
                if ND(tCAone(io,2),2) > (0.5*iy)
                    if iy == 1
                        bCA = [bCA;tCAone(io,:)];
                        bAvar = [bAvar,...
                            (Avar(ismember(CA,tCAone(io,:),'rows'))./2)];
                        brvar = [brvar,...
                            (rvar(ismember(CA,tCAone(io,:),'rows'))./2)];
                    elseif iy == 2
                        tCA = [tCA;tCAone(io,:)];
                        tAvar = [tAvar,...
                            (Avar(ismember(CA,tCAone(io,:),'rows'))./2)];
                        trvar = [trvar,...
                            (rvar(ismember(CA,tCAone(io,:),'rows'))./2)];
                    end
                end
            end
            tCAtwo = CA(ismember(CA(:,2),nodeslower(ib)),:);
            for io = 1:1:size(tCAtwo,1)
                if ND(tCAtwo(io,1),2) > (0.5*iy)
                    if iy == 1
                        bCA = [bCA;tCAtwo(io,:)];
                        bAvar = [bAvar,...
                            (Avar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                        brvar = [brvar,...
                            (rvar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                    elseif iy == 2
                        tCA = [tCA;tCAtwo(io,:)];
                        tAvar = [tAvar,...
                            (Avar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                        trvar = [trvar,...
                            (rvar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                    end
                end
            end
        end
        for ib = 1:1:size(nodesupper,1) % For nodesupper
            tCAone = CA(ismember(CA(:,1),nodesupper(ib)),:);
            for io = 1:1:size(tCAone,1)
                if ND(tCAone(io,2),2) < (0.5*iy)
                    if iy == 1
                        bCA = [bCA;tCAone(io,:)];
                        bAvar = [bAvar,...
                            (Avar(ismember(CA,tCAone(io,:),'rows'))./2)];
                        brvar = [brvar,...
                            (rvar(ismember(CA,tCAone(io,:),'rows'))./2)];
                    elseif iy == 2
                        tCA = [tCA;tCAone(io,:)];
                        tAvar = [tAvar,...
                            (Avar(ismember(CA,tCAone(io,:),'rows'))./2)];
                        trvar = [trvar,...
                            (rvar(ismember(CA,tCAone(io,:),'rows'))./2)];
                    end
                end
            end
            tCAtwo = CA(ismember(CA(:,2),nodesupper(ib)),:);
            for io = 1:1:size(tCAtwo,1)
                if ND(tCAtwo(io,1),2) < (0.5*iy)
                    if iy == 1
                        bCA = [bCA;tCAtwo(io,:)];
                        bAvar = [bAvar,...
                            (Avar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                        brvar = [brvar,...
                            (rvar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                    elseif iy == 2
                        tCA = [tCA;tCAtwo(io,:)];
                        tAvar = [tAvar,...
                            (Avar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                        trvar = [trvar,...
                            (rvar(ismember(CA,tCAtwo(io,:),'rows'))./2)];
                    end
                end
            end
        end
        
    %}
    end
    [tCA,tia,~] = unique(tCA,'rows');
    tAvar = tAvar(tia); trvar = trvar(tia);
    [bCA,bia,~] = unique(bCA,'rows');
    bAvar = bAvar(bia); brvar = brvar(bia);
end

% FUNCTION TO CALCULATE VOLUME FRACTION 
function volFrac = calcVF(NC,CA,rvar,Avar,sel,sidenum,smalldim,xoy,side)
    % Initialize variables
    totalTrussVol = 0;
    
    % Find total volume of members in CA
    for i = 1:size(CA,1)
        % Finding element length from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        % Adding current element to total volume of trusses
        totalTrussVol = totalTrussVol + (L*Avar(i));
    end
    
    % Determine orientation, position of CA
    if xoy == 1
        xdim = smalldim; ydim = sidenum;
        if side == 1
            lone = sel; ltwo = sel/2;
        elseif side == 2
            lone = sel; ltwo = sel;
        end
    elseif xoy == 2
        xdim = sidenum; ydim = smalldim;
        if side == 1
            lone = sel/2; ltwo = sel;
        elseif side == 2
            lone = sel; ltwo = sel;
        end
    end
        
    % Finding average side "thickness" due to differing element radii
    horizrads = [];
    for i = 1:1:size(CA,1)
        for j = 1:1:(xdim-1)
            if ((CA(i,1) + (j*sidenum)) == CA(i,2)) && ...
                    (NC(CA(i,1),2) == lone)
                singl = sel/(sidenum-1);
                x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
                y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
                L = sqrt(((x2-x1)^2)+((y2-y1)^2));
                for k = 1:1:(L/singl)
                    horizrads = [horizrads,rvar(i)];
                end
            elseif ((CA(i,1) - (j*sidenum)) == CA(i,2)) && ...
                    (NC(CA(i,1),2) == lone)
                singl = sel/(sidenum-1);
                x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
                y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
                L = sqrt(((x2-x1)^2)+((y2-y1)^2));
                for k = 1:1:(L/singl)
                    horizrads = [horizrads,rvar(i)];
                end
            end
        end
    end
    vertrads = [];
    for i = 1:1:size(CA,1)
        for j = 1:1:(ydim-1)
            if ((CA(i,1) + j) == CA(i,2)) && (NC(CA(i,1),1) == ltwo)
                singl = sel/(sidenum-1);
                x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
                y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
                L = sqrt(((x2-x1)^2)+((y2-y1)^2));
                for k = 1:1:(L/singl)
                    vertrads = [vertrads,rvar(i)];
                end
            elseif ((CA(i,1) - j) == CA(i,2)) && (NC(CA(i,1),1) == ltwo)
                singl = sel/(sidenum-1);
                x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
                y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
                L = sqrt(((x2-x1)^2)+((y2-y1)^2));
                for k = 1:1:(L/singl)
                    vertrads = [vertrads,rvar(i)];
                end
            end
        end
    end
    thick = mean([mean(horizrads),mean(vertrads)]);
    
    % Calculating volume fraction
    volFrac = totalTrussVol/(2*thick*(sel*(sel/2))); 
end