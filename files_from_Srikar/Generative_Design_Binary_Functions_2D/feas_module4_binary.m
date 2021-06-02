% FOURTH MODULE: The design must contact surrounding unit cells
% contactbool = true for a design that contacts adjacent unit cells, and
% false for a design that does not contact adjacent unit cells
function contactbool = feas_module4_binary(CA,sidenum)
    % Initialize values
    SortedCA = sortrows(CA);   
    contactbool = 0; 

    % Initialize node category vectors
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
end