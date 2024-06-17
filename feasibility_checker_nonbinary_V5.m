function [feasibilityScore] = feasibility_checker_nonbinary_V5(NC,CA,sel,sidenum)
% This function computes the feasibility score for a design 
% Inputs: nodal position matrix NC 
%         Design Connectivity Array CA         
    feasibilityScore = 1;

    % FIRST CONSTRAINT: members only intersect at nodes (no crossing)
    % Sort points from left to right by x-position
    SortedCA = sortrows(CA);
    ND = NC./sel;
    
    % Develop 4xM matrix of line segment endpoint coordinates, where M is 
    %   the number of truss members.  Each row of format (x1,y1,x2,y2),
    %   where point 1 is leftmost, point 2 is rightmost
    PosA = [ND(SortedCA(:,1),1),ND(SortedCA(:,1),2),...
            ND(SortedCA(:,2),1),ND(SortedCA(:,2),2)];
    
    % Loop through each pair of elements
    for i = 1:1:size(PosA,1)
        for j = 1:1:size(PosA,1)
            % Determine whether the given pair of elements intersects
            intersect = findLineSegIntersection([PosA(i,1),PosA(i,2)],...
                        [PosA(i,3),PosA(i,4)],[PosA(j,1),PosA(j,2)],...
                        [PosA(j,3),PosA(j,4)]);
            
            % Throw an error, given an intersection
            if intersect == true
                feasibilityScore = feasibilityScore - 0.1;
                if feasibilityScore < 0.1
                    return
                end
            end
        end
    end
    
    % SECOND CONSTRAINT: Elements (of either the same or different lengths)
    %   cannot overlap
    
    % Loop through each element
    for k = 1:1:size(SortedCA,1)
        % Loop through each element again, to consider each possible pair 
        %   of elements
        for q = 1:1:size(SortedCA,1)
            % Isolate startpoint/endpoint coordinates of both members,
            % calculate their slopes
            A = ND(SortedCA(k,1),:); B = ND(SortedCA(k,2),:);
            C = ND(SortedCA(q,1),:); D = ND(SortedCA(q,2),:);
            mk = (B(2)-A(2))/(B(1)-A(1));
            mq = (D(2)-C(2))/(D(1)-C(1));
            mk = round(mk,4); mq = round(mq,4);
            
            % Check if the same element is being compared twice
            if k == q
                continue
            % Check if the elements' slopes are the same (i.e. they are
            % parallel)
            elseif mk == mq
                % Check if the elements are horizontal and their
                % x-coordinates overlap
                if (mk == 0) && ((C(1) >= A(1)) && (C(1) < B(1)))
                    % Check if the elements' y-coordinates overlap
                    if C(2) == A(2)
                        feasibilityScore = feasibilityScore - 0.1;
                        if feasibilityScore < 0.1
                            return
                        end
                    end
                % Check if the elements are vertical and their coordinates
                % overlap
                elseif isinf(mk) && ((C(2) >= A(2)) && (C(2) < B(2)))
                    % Check if the elements' x-coordinates overlap
                    if C(1) == A(1)
                        feasibilityScore = feasibilityScore - 0.1;
                        if feasibilityScore < 0.1
                            return
                        end
                    end
                else
                    if k < q
                        t1 = (C(1)-A(1))/(B(1)-A(1));
                        t2 = (C(2)-A(2))/(B(2)-A(2));
                        t3 = (D(1)-A(1))/(B(1)-A(1));
                    elseif k > q
                        t1 = (A(1)-C(1))/(D(1)-C(1));
                        t2 = (A(2)-C(2))/(D(2)-C(2));
                        t3 = (B(1)-C(1))/(D(1)-C(1));
                    end
                    % Check if the diagonal elements overlap
                    if (t1 == t2)
                        if (t1 >= 0) && (t1 < 1)
                            feasibilityScore = feasibilityScore - 0.1;
                            if feasibilityScore < 0.1
                                return
                            end
                        elseif (t3 >= 0) && (t3 < 1)
                            feasibilityScore = feasibilityScore - 0.1;
                            if feasibilityScore < 0.1
                                return
                            end
                        end
                    end
                end
            end
        end
    end
    
    % THIRD CONSTRAINT: The design must have intracell connectivity
    flowboolassn = 0; flowbool = 0;   
    
    % Find unique startpoints and endpoints
    startpoints = unique(SortedCA(:,1)); endpoints = unique(SortedCA(:,2));
    startunique = setdiff(startpoints,endpoints);
    endunique = setdiff(endpoints,startpoints);
    
    % Check if there are any unique start- or endpoints at all
    if isempty(startunique) || isempty(endunique)
        flowbool = 1;
    else
        % For each unique startpoint, determine if connectivity exists to
        % all unique endpoints
        endreached = zeros(1,length(endunique));
        for v = 1:1:length(startunique)
            mCA = SortedCA(SortedCA(:,1)==startunique(v),:);
            while size(mCA,1) < size(SortedCA,1)
                nmCA = mCA;
                noneleft = 1;
                for z = 1:1:size(mCA,1)
                    tCA = SortedCA(SortedCA(:,1)== (mCA(z,2)),:);
                    repeatrows = [];
                    for q = 1:1:size(tCA,1)
                        if ismember(tCA(q,:),mCA,'rows')
                            repeatrows = [repeatrows,q];
                        end
                    end
                    tCA(repeatrows,:) = [];
                    if ~isempty(tCA)
                        noneleft = 0;
                    end
                    nmCA = [nmCA;tCA];
                end
                nmCA = unique(nmCA,'rows');
                flowdet = ismember(endunique,nmCA(:,2));
                % If the current unique startpoint leads to all unique
                % endpoints, the design is intraconnected
                if all(flowdet) == true
                    flowbool = 1;
                    break
                % If this unique startpoint only connects to one unique
                % endpoint, this portion of the design could be isolated
                % from the rest, depending on whether any other flowpaths
                % terminate in the current one
                elseif (noneleft == 1) && (sum(flowdet) == 1)
                    usednodes = unique(nmCA);
                    extflowdet = ismember(usednodes,SortedCA(:,2));
                    if ~any(extflowdet)
                        flowbool = 0;
                        disp('here');
                        break
                    else
                        endreached(flowdet==1) = endreached(flowdet==1)+1;
                        break
                    end
                % If there are no longer any members in the possible
                % flowpaths, this unique startpoint doesn't connect to all
                % unique endpoints  
                elseif noneleft == 1
                    endreached(flowdet==1) = endreached(flowdet==1) + 1;
                    break
                end
                mCA = nmCA;
            end
        end
        %{
        % If a single startpoint doesn't lead to all endpoints, determine
        % if all startpoints cumulatively lead to all endpoints
        if all(endreached) && any(endreached>1)
            flowbool = 1;
        end
        %}
    end
    
    % Score design based on number of unconnected nodes
    if flowbool ~= 1
        penalty = max(length(startunique),length(endunique));
        feasibilityScore = feasibilityScore - (0.1*penalty);
        if feasibilityScore < 0.1
            return
        end
    end
    
    % FOURTH CONSTRAINT: The design must have intercell connectivity
    % Initialize node category vectors
    contactboolx = 0; contactbooly = 0;
    usedpoints = unique(SortedCA);
    lenodes = 2:1:(sidenum-1); 
    renodes = ((sidenum^2)-sidenum+2):1:((sidenum^2)-1);
    tenodes = (2*sidenum):sidenum:((sidenum^2)-sidenum);
    benodes = (sidenum+1):sidenum:((sidenum^2)-(2*sidenum)+1);
    cornernodes = [1,sidenum,((sidenum^2)-sidenum+1),(sidenum^2)];
    lrnodepairs = [lenodes;renodes]';
    tbnodepairs = [benodes;tenodes]';
    
    % Determine if design has at least one repeated external contact point
    for y = 1:1:length(usedpoints)
        otherpoints = usedpoints([1:(y-1),(y+1):end]);
        if ismember(usedpoints(y),cornernodes)
            ocnode = usedpoints(y);
            if ocnode == 1
                lrcnodes = [((sidenum^2)-sidenum+1),(sidenum^2)];
                tbcnodes = [sidenum,(sidenum^2)];
            elseif ocnode == sidenum
                lrcnodes = [((sidenum^2)-sidenum+1),(sidenum^2)];
                tbcnodes = [1,((sidenum^2)-sidenum+1)];
            elseif ocnode == ((sidenum^2)-sidenum+1)
                lrcnodes = [1,sidenum];
                tbcnodes = [sidenum,(sidenum^2)];
            elseif ocnode == (sidenum^2)
                lrcnodes = [1,sidenum];
                tbcnodes = [1,((sidenum^2)-sidenum+1)];
            end
            lrcornercontact = ismember(lrcnodes,otherpoints);
            tbcornercontact = ismember(tbcnodes,otherpoints);
            if any(lrcornercontact) && any(tbcornercontact)
                contactboolx = 1; contactbooly = 1;
            elseif any(lrcornercontact)
                contactboolx = 1;
            elseif any(tbcornercontact)
                contactbooly = 1;
            end
        elseif ismember(usedpoints(y),lrnodepairs)
            [i,j] = find(lrnodepairs == usedpoints(y));
            if j == 1
                if ismember(lrnodepairs(i,2),otherpoints)
                    contactboolx = 1;
                end
            elseif j == 2
                if ismember(lrnodepairs(i,1),otherpoints)
                    contactboolx = 1;
                end
            end
        elseif ismember(usedpoints(y),tbnodepairs)
            [i,j] = find(tbnodepairs == usedpoints(y));
            if j == 1
                if ismember(tbnodepairs(i,2),otherpoints)
                    contactbooly = 1;
                end
            elseif j == 2
                if ismember(tbnodepairs(i,1),otherpoints)
                    contactbooly = 1;
                end
            end    
        end
    end
    
    % Determine if intercell connectivity exists in both x and y
    if (contactboolx == 1) && (contactbooly == 1)
        contactbool = 1;
    else
        contactbool = 0;
    end
    
    % Score design unilaterally based on absence of contact points
    if contactbool ~= 1
        feasibilityScore = feasibilityScore - 0.5;
        if feasibilityScore < 0.1
            return
        end
    end
end

% FUNCTION TO DETERMINE PRESENCE OF INTERSECTION
function intersect = findLineSegIntersection(p1,q1,p2,q2)
    if (findOrientation(p1,q1,p2) ~= findOrientation(p1,q1,q2))&&...
            (findOrientation(p2,q2,p1) ~= findOrientation(p2,q2,q1)) 
        if isequal(p1,p2) || isequal(q1,q2) || ...
                isequal(p1,q2) || isequal(q1,p2)
            % Intersection due to shared start or end point
            intersect = false;
        else
            intersect = true;
        end
    else
        % Intersection not present
        intersect = false;
    end
end

% FUNCTION TO CALCULATE ORIENTATION FROM 3 POINTS
function orientation = findOrientation(p,q,r)
    val = ((q(2)-p(2))*(r(1)-q(1)))-((q(1)-p(1))*(r(2)-q(2)));
    if (abs(val) < (10^-5)) && (val > 0)
        orientation = 0;
    elseif val > (10^-5)
        orientation = 1;
    else
        orientation = 2;
    end
end

