% This function calculates the instances of intersection in a design, and
% outputs a vector with three columns.  The first two columns represent the
% connectivity of a member (the two nodes it connects to).  The third
% column represents the number of times that member intersects with other
% members.  The number of rows in the output signifies the numebr of
% intersecting members present in the design.
function outputMemLog = intersectLogger(NC,CA_des)
    % Initialize values
    SortedCA = sortrows(CA_des);
    
    % Develop 4xM matrix of line segment endpoint coordinates, where M is 
    %   the number of truss members.  Each row of format (x1,y1,x2,y2),
    %   where point 1 is leftmost, point 2 is rightmost
    PosA = [NC(SortedCA(:,1),1),NC(SortedCA(:,1),2),...
            NC(SortedCA(:,2),1),NC(SortedCA(:,2),2)];
    
    % Loop through each pair of elements to find intersecting pairs
    memberLog = [];
    for i = 1:1:size(PosA,1)
        for j = 1:1:size(PosA,1)
            % Determine whether the given pair of elements intersects
            intersect = findLineSegIntersection([PosA(i,1),PosA(i,2)],...
                        [PosA(i,3),PosA(i,4)],[PosA(j,1),PosA(j,2)],...
                        [PosA(j,3),PosA(j,4)]);
            
            % Given an intersection, log the two intersecting members
            if intersect == true
                memberLog = [memberLog;SortedCA(i,:);SortedCA(j,:)];
            end
        end
    end
    
    % Find number of times each member intersects another
    uniquemems = unique(memberLog,'rows');
    totalcounts = zeros(size(uniquemems,1),1);
    for q = 1:1:size(uniquemems,1)
        intmem = uniquemems(q);
        % ismember for each member
        presencelog = ismember(memberLog,uniquemems(q,:),'rows');
        totalcounts(q) = sum(presencelog);
    end
    totalcounts = totalcounts./2;
    outputMemLog = [uniquemems,totalcounts];
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
