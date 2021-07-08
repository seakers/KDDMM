% Symmetry Heuristic V3
% This function determines whether a design is fully symmetric (behaves the
% same when rotated 90 degrees four times), flip symmetric (behaves the
% same when rotated 180 degrees twice), or asymmmetric.  The score assigned
% to full symmetry is 2, the score assigned to flip symmetry is 1, and teh
% score assigned to asymmetry is 0
function symmScore = symmHeuristic_2D_V3(CA,sel,sidenum)
    % Initialize score
    symmScore = 2;

    % Generate vector with nodal coordinates
    NC = generateNC(sel,sidenum);
    
    % Check for full symmetry
    for i = 1:1:size(CA,1)
        % Find the rotation of both endpoints by 90 degrees (shifted back
        % to the first quadrant)
        ry1 = round((NC(CA(i,1),1)),6); 
        ry2 = round((NC(CA(i,2),1)),6);
        rx1 = round((-NC(CA(i,1),2) + sel),6); 
        rx2 = round((-NC(CA(i,2),2) + sel),6);
        
        % Find nodal indices of reflected points
        Ax = find(NC(:,1)==rx1); Ay = find(NC(:,2)==ry1);
        A = intersect(Ax,Ay);
        Bx = find(NC(:,1)==rx2); By = find(NC(:,2)==ry2);
        B = intersect(Bx,By);
        
        % Determine whether reflected members exist in CA (and demerit
        % design if not)
        if (ismember([A,B],CA,'rows') == false) && ...
           (ismember([B,A],CA,'rows') == false)
            symmScore = 0;
            break;
        end
    end
    
    % Check for flip symmetry
    if symmScore == 0
        symmScore = 1; fsChecker = 2;
        
        % Vertical flip symmetry
        for i = 1:1:size(CA,1)
            % Find the flip of both endpoints about the x axis (shifted 
            % back to the first quadrant)
            rx1 = round((NC(CA(i,1),1)),6); 
            rx2 = round((NC(CA(i,2),1)),6);
            ry1 = round((-NC(CA(i,1),2) + sel),6); 
            ry2 = round((-NC(CA(i,2),2) + sel),6);

            % Find nodal indices of x-reflected points
            Fx = find(NC(:,1)==rx1); Fy = find(NC(:,2)==ry1);
            F = intersect(Fx,Fy);
            Gx = find(NC(:,1)==rx2); Gy = find(NC(:,2)==ry2);
            G = intersect(Gx,Gy);

            % Determine whether reflected members exist in CA (and demerit
            % design if not)
            if (ismember([F,G],CA,'rows') == false) && ...
               (ismember([G,F],CA,'rows') == false)
                fsChecker = fsChecker - 1;
                break;
            end
        end
        
        % Horizontal flip symmetry
        for i = 1:1:size(CA,1)
            % Find the flip of both endpoints about the x axis (shifted 
            % back to the first quadrant)
            ry1 = round((NC(CA(i,1),1)),6); 
            ry2 = round((NC(CA(i,2),1)),6);
            rx1 = round((-NC(CA(i,1),2) + sel),6); 
            rx2 = round((-NC(CA(i,2),2) + sel),6);

            % Find nodal indices of y-reflected points
            Jx = find(NC(:,1)==rx1); Jy = find(NC(:,2)==ry1);
            J = intersect(Jx,Jy);
            Dx = find(NC(:,1)==rx2); Dy = find(NC(:,2)==ry2);
            D = intersect(Dx,Dy);

            % Determine whether reflected members exist in CA (and demerit
            % design if not)
            if (ismember([J,D],CA,'rows') == false) && ...
               (ismember([D,J],CA,'rows') == false)
                fsChecker = fsChecker - 1;
                break;
            end
        end
        
        % Check for either flip symmetry
        if fsChecker == 0
            symmScore = 0;
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