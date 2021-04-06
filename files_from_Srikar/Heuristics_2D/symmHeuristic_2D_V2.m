% Symmetry Heuristic
% This function determines whether a design is symmetric about the y=x
% line, and scores it based on the number of members that are not symmetric
function symmScore = symmHeuristic_2D_V2(CA,sel,sidenum)
    % Initialize scores
    symmScore = 1;

    % Generate vector with nodal coordinates
    NC = generateNC(sel,sidenum);
    
    % Loop through each member
    for i = 1:1:size(CA,1)
        % Find the rotation of both endpoints by 90 degrees (shifted back
        % to the first quadrant)
        ry1 = NC(CA(i,1),1); ry2 = NC(CA(i,2),1);
        rx1 = -NC(CA(i,1),2) + sel; rx2 = -NC(CA(i,2),2) + sel;
        
        % Find nodal indices of reflected points
        Ax = find(NC(:,1)==rx1); Ay = find(NC(:,2)==ry1);
        A = intersect(Ax,Ay);
        Bx = find(NC(:,1)==rx2); By = find(NC(:,2)==ry2);
        B = intersect(Bx,By);
        
        % Determine whether reflected members exist in CA (and demerit
        % design if not)
        if (ismember([A,B],CA,'rows') == false) && ...
           (ismember([B,A],CA,'rows') == false)
            symmScore = symmScore-0.1;
            if symmScore < 0.1
                return
            end
        end
    end
end

%-------%
% FUNCTION TO GENERATE NODAL COORDINATES BASED ON GRID SIZE
function NC = generateNC(sel,sidenum)
    notchvec = linspace(0,1,sidenum);
    NC = [];
    for i = 1:1:sidenum
        for j = 1:1:sidenum
            NC = [NC;notchvec(i),notchvec(j)];
        end
    end
    NC = sel.*NC;
end