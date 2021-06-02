% INTERCONNECTIVITY CHECKER
% This function explicitly checks if unit cell designs are conencted within 
% each other as well as connected to adjacent unit cells, and returns a
% binary result
function icBool = interconnectivityChecker(sidenum,CA)
    % Initialize values
    SortedCA = sortrows(CA);

    % FIRST MODULE: The entire design must be interconnected
    % Find unique startpoints and endpoints
    startpoints = unique(SortedCA(:,1)); endpoints = unique(SortedCA(:,2));
    startunique = setdiff(startpoints,endpoints);
    endunique = setdiff(endpoints,startpoints);
    
    % For each unique startpoint, determine if connectivity exists to at
    % least one unique endpoint
    flowbool = 0;
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
    
    % SECOND MODULE: The design must contact surrounding unit cells
    % Initialize node category vectors
    contactbool = 0;
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
                contactbool = 1;
            end
        elseif ismember(usedpoints(y),lrnodepairs)
            [i,j] = find(lrnodepairs == usedpoints(y));
            if j == 1
                if ismember(lrnodepairs(i,2),otherpoints)
                    contactbool = 1;
                end
            elseif j == 2
                if ismember(lrnodepairs(i,1),otherpoints)
                    contactbool = 1;
                end
            end
        elseif ismember(usedpoints(y),tbnodepairs)
            [i,j] = find(tbnodepairs == usedpoints(y));
            if j == 1
                if ismember(tbnodepairs(i,2),otherpoints)
                    contactbool = 1;
                end
            elseif j == 2
                if ismember(tbnodepairs(i,1),otherpoints)
                    contactbool = 1;
                end
            end    
        end
    end
    
    % Score design unilaterally based on absence of contact points
    if contactbool ~= 1
        feasibilityScore = feasibilityScore - 0.5;
        if feasibilityScore < 0.1
            return
        end
    end
    
    % DETERMINE OUTPUT BOOLEAN
    if (contactbool == 1) && (flowbool == 1)
        icBool = 1;
    else
        icBool = 0;
    end
end