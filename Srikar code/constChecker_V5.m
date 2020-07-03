% Constraint Checker
% This function is intended to check that a given truss design is
% legitimate (not necessarily stable, but legitimate within the standard
% nodal framework and related assumptions)
function [constraintScore,constraintBool] = constChecker_V5(NC,CA)
    constraintScore = 1;

    % FIRST CONSTRAINT: members only intersect at nodes (no crossing)
    % (source: https://www.geeksforgeeks.org/given-a-set-of-line-segments
    % -find-if-any-two-segments-intersect/)
    
    % Sort points from left to right by x-position
    SortedCA = sortrows(CA);
    
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
                constraintScore = constraintScore - 0.1;
                if constraintScore < 0.1
                    return
                end
                %{
                D = ['The element from (',num2str(PosA(i,1)),',',...
                     num2str(PosA(i,2)),') to (',num2str(PosA(i,3)),...
                     ',',num2str(PosA(i,4)),') intersects with the'...
                     ' element from (',num2str(PosA(j,1)),',',...
                     num2str(PosA(j,2)),') to (',num2str(PosA(j,3)),',',...
                     num2str(PosA(j,4)),')'];
                disp(D);
                %}
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
            % Check if both elements share a common startpoint
            if (NC(SortedCA(k,1),1) == NC(SortedCA(q,1),1)) && ...
                (NC(SortedCA(k,1),2) == NC(SortedCA(q,1),2))
                % Check if both elements have the same slope (and reject 
                %    the design if so)
                mk = (NC(SortedCA(k,2),2)-NC(SortedCA(k,1),2))/...
                     (NC(SortedCA(k,2),1)-NC(SortedCA(k,1),1));
                mq = (NC(SortedCA(q,2),2)-NC(SortedCA(q,1),2))/...
                     (NC(SortedCA(q,2),1)-NC(SortedCA(q,1),1));
                if mk == mq
                   constraintScore = constraintScore - 0.1;
                   if constraintScore < 0.1
                       return
                   end
                   %{
                   D = ['The element from (',num2str(PosA(k,1)),',',...
                     num2str(PosA(k,2)),') to (',num2str(PosA(k,3)),...
                     ',',num2str(PosA(k,4)),') overlaps with the'...
                     ' element from (',num2str(PosA(q,1)),',',...
                     num2str(PosA(q,2)),') to (',num2str(PosA(q,3)),',',...
                     num2str(PosA(q,4)),')'];
                   disp(D);
                   %}
                end
            end
        end
    end
    
    % Assign value to constraint boolean
    constraintBool = true;
    if constraintScore < 1
        constraintBool = false;
    end
end

% FUNCTION TO DETERMINE PRESENCE OF INTERSECTION (FOR CONSTRAINT #1)
% (source: https://www.geeksforgeeks.org/check-if-two-given-line-segments
% -intersect/)
% This boolean function determines whether two line segments intersect,
% given their endpoints as inputs
function intersect = findLineSegIntersection(p1,q1,p2,q2)
    if (findOrientation(p1,q1,p2) ~= findOrientation(p1,q1,q2))&&...
            (findOrientation(p2,q2,p1) ~= findOrientation(p2,q2,q1)) 
        if isequal(p1,p2) || isequal(q1,q2) || ...
                isequal(p1,q2) || isequal(q1,p2)
            intersect = false; 
        else
            intersect = true;
        end
    else
        intersect = false;
    end
end

% FUNCTION TO CALCULATE ORIENTATION FROM 3 POINTS (FOR CONSTRAINT #1)
% (source: https://www.geeksforgeeks.org/orientation-3-ordered-points/)
% This function finds the orientation of an ordered triplet (p, q, r)
% The function returns one of three following values:
% 0 --> p, q and r are colinear 
% 1 --> Clockwise 
% 2 --> Counterclockwise 
function orientation = findOrientation(p,q,r)
    val = ((q(2)-p(2))*(r(1)-q(1)))-((q(1)-p(1))*(r(2)-q(2)));
    if val == 0
        orientation = 0;
    elseif val > 0
        orientation = 1;
    else
        orientation = 2;
    end
end

