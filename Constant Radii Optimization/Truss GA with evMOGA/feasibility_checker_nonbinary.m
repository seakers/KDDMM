function [feasibilityScore] = feasibility_checker_nonbinary(NC,CA_des)
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
                % If the same element is being compared twice, move on
                if k == q
                    continue
                elseif mk == mq
                   feasibilityScore = feasibilityScore - 0.1;
                   if feasibilityScore < 0.1
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

end

