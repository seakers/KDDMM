% THIRD MODULE: The entire design must be connected within itself
% (intracell connecvitiy)
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
    
    % Check if there are any unique start- or endpoints at all
    if isempty(startunique) || isempty(endunique)
        flowbool = 1;
        return
    else
        % Determine whether there are startpoints or endpoints that are
        % only used once in the CA
        usednodes = unique(SortedCA);
        freqvec = [usednodes,histc(SortedCA(:),usednodes)];
        singlestarts = []; singleends = [];
        for i = 1:1:length(startunique)
            if freqvec((freqvec(:,1) == startunique(i)),2) == 1
                singlestarts = [singlestarts,startunique(i)];
            end
        end
        for j = 1:1:length(endunique)
            if freqvec((freqvec(:,1) == endunique(j)),2) == 1
                singleends = [singleends,endunique(j)];
            end
        end
        if ~isempty(singlestarts) && ~isempty(singleends)
            for k = 1:1:length(singlestarts)
                for m = 1:1:length(singleends)
                    member = [singlestarts(k),singleends(m)];
                    if ismember(member,SortedCA,'rows')
                        flowbool = 0;
                        return
                    end
                end
            end
        end
        
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
                nmCA = unique(nmCA,'rows');
                flowdet = ismember(endunique,nmCA(:,2));
                if any(flowdet) == true
                    flowbool = 1;
                    break
                end
                mCA = nmCA;
            end
            if flowbool == 1
                break
            end
        end
    end
end