% THIRD MODULE: The entire design must be interconnected
% flowbool = true for a unit cell that's fully connected within itself, and
% false for a unit cell that is not
function flowbool = feas_module3_binary(CA)
    % Initialize values
    SortedCA = sortrows(CA);
    flowbool = 0;

    % Find unique startpoints and endpoints
    startpoints = unique(SortedCA(:,1)); endpoints = unique(SortedCA(:,2));
    startunique = setdiff(startpoints,endpoints);
    endunique = setdiff(endpoints,startpoints);
    
    % For each unique startpoint, determine if connectivity exists to at
    % least one unique endpoint
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