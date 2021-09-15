% FUNCTION TO GENERATE 5x5 ROTATIONALLY SYMMETRIC DESIGNS FROM 3x3 DESIGNS
function sCABasket_5x5 = symm5x5From3x3(CABasket_3x3,NC5,NC3,sel,sidenum)
    % Generate symmetric 5x5 design from each 3x3 design
    sCABasket_5x5 = {};
    for i = 1:1:length(CABasket_3x3)
        CA_3x3 = cell2mat(CABasket_3x3(i));
        CA_5x5 = [];
        
        % Loop through each member in the 3x3 design
        for j = 1:1:size(CA_3x3,1)
            % Input stock design (Pos. 1) to CA_5x5 using NC5 coordinates
            x1 = NC3(CA_3x3(j,1),1); x2 = NC3(CA_3x3(j,2),1);
            y1 = NC3(CA_3x3(j,1),2); y2 = NC3(CA_3x3(j,2),2);
            Ax = find(NC5(:,1)==x1); Ay = find(NC5(:,2)==y1);
            A = intersect(Ax,Ay);
            Bx = find(NC5(:,1)==x2); By = find(NC5(:,2)==y2);
            B = intersect(Bx,By);
            CA_5x5 = [CA_5x5;[A,B]];
            
            % Rotate each member three times, input rotated member to
            % CA_5x5 using NC5 location
            for k = 1:1:3
                ry1 = round(4*x1,2)/4; ry2 = round(4*x2,2)/4; 
                rx1 = round(4*(-y1 + (sel)),2)/4; 
                rx2 = round(4*(-y2 + (sel)),2)/4;
                Qx = find(NC5(:,1)==rx1); Qy = find(NC5(:,2)==ry1);
                Q = intersect(Qx,Qy);
                Dx = find(NC5(:,1)==rx2); Dy = find(NC5(:,2)==ry2);
                D = intersect(Dx,Dy);
                CA_5x5 = [CA_5x5;[Q,D]];
                x1 = rx1; x2 = rx2; y1 = ry1; y2 = ry2;
            end
        end
        
        % Filter out repeated members, sort CA
        sortCA = sortrows((sort(CA_5x5'))');
        sfCA = unique(sortCA,'rows');
        
        % Check for repeatability (and amend design with full frame if not)
        repBool = repChecker(sfCA,5);
        if repBool ~= 1
            leftCA = [(1:1:(sidenum-1))',(2:1:sidenum)'];
            topCA = [(sidenum:sidenum:((sidenum^2)-sidenum))',...
                      ((2*sidenum):sidenum:(sidenum^2))'];
            rightCA = [(((sidenum^2)-sidenum+1):1:((sidenum^2)-1))',...
                       (((sidenum^2)-sidenum+2):1:(sidenum^2))'];    
            bottomCA = [(1:sidenum:((sidenum^2)-(2*sidenum)+1))',...
                        ((sidenum+1):sidenum:((sidenum^2)-sidenum+1))'];
            sfCA = [sfCA;leftCA;topCA;rightCA;bottomCA];
            sfCA = unique(sfCA,'rows');
        end
        
        % Add design to collection
        sCABasket_5x5(i) = {sfCA};
    end
end