function [feasibilityScore] = feasibility_checker_nonbinary_V4(NC,CA_des)
% This function computes the feasibility score for a design 
% Inputs: nodal position matrix NC 
%         Design Connectivity Array CA_des         
    feasibilityScore = 1;

    % FIRST CONSTRAINT: members only intersect at nodes (no crossing)
    % Sort points from left to right by x-position
    SortedCA = sortrows(CA_des);
    
    % Develop 4xM matrix of line segment endpoint coordinates, where M is 
    %   the number of truss members.  Each row of format (x1,y1,x2,y2),
    %   where point 1 is leftmost, point 2 is rightmost
    PosA = [NC(SortedCA(:,1),1),NC(SortedCA(:,1),2),...
            NC(SortedCA(:,2),1),NC(SortedCA(:,2),2)];
    
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
            A = NC(SortedCA(k,1),:); B = NC(SortedCA(k,2),:);
            C = NC(SortedCA(q,1),:); D = NC(SortedCA(q,2),:);
            mk = (B(2)-A(2))/(B(1)-A(1));
            mq = (D(2)-C(2))/(D(1)-C(1));
            
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
                    t1 = (C(1)-A(1))/(B(1)-A(1));
                    t2 = (C(2)-A(2))/(B(2)-A(2));
                    t3 = (D(1)-A(1))/(B(1)-A(1));
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
    flowbool = 0;    
    
    % Find unique startpoints and endpoints
    startpoints = unique(SortedCA(:,1)); endpoints = unique(SortedCA(:,2));
    startunique = setdiff(startpoints,endpoints);
    endunique = setdiff(endpoints,startpoints);
    
    if isempty(startunique) || isempty(endunique)
        flowbool = 1;
    else
        % For each unique startpoint, determine if connectivity exists to
        % at least one unique endpoint
        for v = 1:1:length(startunique)
            mCA = SortedCA(SortedCA(:,1)==startunique(v),:);
            while size(mCA,1) < size(SortedCA,1)
                nmCA = mCA;
                for z = 1:1:size(mCA,1)
                    tCA = SortedCA(SortedCA(:,1)== (mCA(z,2)),:);
                    nmCA = [nmCA;tCA];
                end
                for x = 1:1:size(nmCA,1)
                    flowdet = ismember(endunique,nmCA(:,2));
                    if any(flowdet)
                        flowbool = 1;
                        break
                    end
                end
                if flowbool == 1
                    break
                end
                mCA = nmCA;
            end
            if flowbool == 1
                break
            end
        end
    end
    
    % Score design based on number of unconnected nodes
    if flowbool ~= 1
        penalty = max(length(startunique,endunique));
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
            cnodes = cornernodes(cornernodes~=usedpoints(y));
            cornercontact = ismember(cnodes,otherpoints);
            if any(cornercontact)
                contactboolx = 1; contactbooly = 1;
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

