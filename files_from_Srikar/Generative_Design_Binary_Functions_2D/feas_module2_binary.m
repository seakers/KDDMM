% SECOND CONSTRAINT: Elements (of either the same or different lengths)
% nooverlapbool = true for a design with no overlapping members and false
% for a design that contains at least one instance of overlapping members
function nooverlapbool = feas_module2_binary(CA,NC)
    % Initialize values
    SortedCA = sortrows(CA);
    nooverlapbool = 1;

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
                   nooverlapbool = 0;
                   return
                end
            % Check if both elements share a common endpoint    
            elseif (NC(SortedCA(k,2),1) == NC(SortedCA(q,2),1)) && ...
                   (NC(SortedCA(k,2),2) == NC(SortedCA(q,2),2))    
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
                   nooverlapbool = 0;
                   return
                end
            end
        end
    end
end