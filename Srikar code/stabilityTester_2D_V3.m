% FUNCTION TO TEST 2D TRUSS STABILITY 
function [stabilityBool,stabilityScore] = ...
                                       stabilityTester_2D_V3(sidenum,CA,NC)
    % Initialize stability score
    stabilityScore = 1;

    % Add up counters based on nodal connectivities
    [N,~] = histcounts(CA,size(NC,1));
    
    % First stability check: number of "holes" (unconnected nodes) in truss
    %   should be less than or equal to [(number of side nodes) - 2]
    zeros = find(~N);
    if length(zeros) > (sidenum-2)
        stabilityScore = stabilityScore - 0.1;
        if stabilityScore < 0.1
            return
        end
    end
    
    % Second stability check: nodes with connections are connected to at
    %   least three other nodes apiece (except for the corner nodes)
    Ns = N([2:(sidenum-1),(sidenum+1):((sidenum^2)-sidenum),...
         ((sidenum^2)-(sidenum-2)):(sidenum^2)-1]);
    Nnz = Ns(Ns>0);
    for a = 1:1:length(Nnz)
       if (Nnz(a) == 1) || (Nnz(a) == 2)
           stabilityScore = stabilityScore - 0.1;
           if constraintScore < 0.1
               return
           end
       end
    end
    
    % Third stability check: corner nodes have at least two connections
    Nc = N([1,sidenum,((sidenum^2)-(sidenum-1)),sidenum^2]);
    for a = 1:1:length(Nc)
       if Nc(a) == 1
           stabilityScore = stabilityScore - 0.2;
           if stabilityScore < 0.1
               return
           end
       end
    end
    
    % Fourth stability check: at least one diagonal member present
    nodiags = true;
    for i = 1:1:size(CA,1)
        if CA(i,1)+sidenum == CA(i,2)
            nodiags = true;
        elseif CA(i,1)-sidenum == CA(i,2)
            nodiags = true;
        elseif CA(i,1)+1 == CA(i,2)
            nodiags = true;
        elseif CA(i,1)-1 == CA(i,2)
            nodiags = true;
        else
            nodiags = false;
            break
        end
    end
    if nodiags == true
        stabilityScore = stabilityScore - 0.2;
        if stabilityScore < 0.1
            return
        end
    end
    
    % Assign value to stability boolean
    stabilityBool = true;
    if stabilityScore < 1
        stabilityBool = false;
    end
end