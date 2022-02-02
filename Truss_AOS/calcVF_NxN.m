% FUNCTION TO CALCULATE VOLUME FRACTION 
function volFrac = calcVF_NxN(CA,r,sel,sidenum)
    rvar = r.*ones(1,size(CA,1));

    % Generate vector with nodal coordinates
    NC = generateNC(sel,sidenum);
    
    % Calculate Avar & modify for edge members
    Avar = pi.*(rvar.^2); % Cross-sectional areas of truss members
    Avar = modifyAreas(Avar,CA,NC,sidenum);

    % Calculate cumulative volume of all individual members
    totalTrussVol = 0;
    totl = sel;
    for i = 1:size(CA,1)
        % Finding element length from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        % Adding current element to total volume of members
        totalTrussVol = totalTrussVol + (L*Avar(i));
    end
    
    % Modify total member volume based on intersecting and overlapping
    % members
    totalTrussVol = modVol(NC,CA,totalTrussVol,r);
    %disp(totalTrussVol);
    
    % Modify total member volume based on overlaps at nodes
    totalTrussVol = subNodOLVol(NC,CA,totalTrussVol,rvar);
    
    % Finding average side "thickness" due to differing element radii
    thick = mean(rvar);
    
    % Calculating volume fraction (using a solid square with 2*(avg 
    %   thickness) as a baseline)
    volFrac = totalTrussVol/(2*thick*(totl^2)); 
end

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

% FUNCTION TO MODIFY AREAS FOR EDGE MEMBERS
function Avar = modifyAreas(Avar,CA,NC,sidenum)
    % Identify edge nodes
    edgenodes = [(1:1:sidenum),((2*sidenum):sidenum:...
                 ((sidenum^2)-sidenum)),((sidenum+1):sidenum:...
                 ((sidenum^2)-(2*sidenum)+1)),(((sidenum^2)-sidenum+1):...
                 1:(sidenum^2))];
             
    % Identify members connecting solely to edge nodes
    edgeconn1 = ismember(CA(:,1),edgenodes);
    edgeconn2 = ismember(CA(:,2),edgenodes);
    edgeconnectors = (edgeconn1 & edgeconn2);
    
    % Isolate edge members based on angle
    edgelogical = [edgeconnectors,edgeconnectors];
    CAedgenodes = CA.*(edgelogical);
    CAedgenodes = CAedgenodes(any(CAedgenodes,2),:);
    x1 = NC(CAedgenodes(:,1),1); x2 = NC(CAedgenodes(:,2),1);
    y1 = NC(CAedgenodes(:,1),2); y2 = NC(CAedgenodes(:,2),2);
    L = sqrt(((x2-x1).^2)+((y2-y1).^2));
    angles = rad2deg(abs(acos((x2-x1)./L)));
    CAedgy = [];
    for i = 1:1:size(CAedgenodes,1)
        if (angles(i) == 0) || (angles(i) == 90)
            CAedgy = [CAedgy;CAedgenodes(i,:)];
        end
    end
    
    % Find and modify areas belonging to edge members
    if isempty(CAedgy)
        % Do nothing
    else
        edgemembers = ismember(CA,CAedgy,'rows');
        selectAreas = Avar'.*edgemembers;
        k = find(selectAreas);
        Avar(k) = Avar(k)./2;
    end
end

% FUNCTION TO MODIFY TOTAL TRUSS VOLUME BASED ON INTERSECTIONS
function tTV = modVol(NC,CA,tTV,r)
    %disp(tTV);
    % NOTE: Both the crossing and overlap modules are not capable of
    % explicitly ignoring an instance after it has been identified the 
    % first time. The overlap volumes subtracted herein are half of the 
    % actual total overlap volume for each instance; for a single crossing 
    % the statements will be run twice
    
    % Subtract duplicate volumes of crossing elements
    SortedCA = sortrows(CA);
    
    % Develop 4xM matrix of line segment endpoint coordinates, where M is 
    %   the number of truss members.  Each row of format (x1,y1,x2,y2),
    %   where point 1 is leftmost, point 2 is rightmost
    PosA = [NC(SortedCA(:,1),1),NC(SortedCA(:,1),2),...
            NC(SortedCA(:,2),1),NC(SortedCA(:,2),2)];
    %{
    % Loop through each pair of elements
    for i = 1:1:size(PosA,1)
        for j = 1:1:size(PosA,1)
            % Determine whether the given pair of elements intersects
            [intersect,alpha] = findLineSegIntersection([PosA(i,1),...
                        PosA(i,2)],[PosA(i,3),PosA(i,4)],[PosA(j,1),...
                        PosA(j,2)],[PosA(j,3),PosA(j,4)]);
            
            % Subtract overlap volume, given an intersection
            if (intersect == true) && (alpha ~= 0)
                %tTV = tTV-((8*(r^3))/(3*abs(sin(alpha))));
                tTV = tTV - ((2*pi*(r^3))/3);
                disp(tTV);
                fprintf('[i,j] = [%.5f, %.5f]\n',[i,j])
            end
        end
    end
    %}
    
    % Subtract duplicate volumes of overlapping elements
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
                % If the same element is being compared to itself, move on
                if k == q
                    continue
                elseif mk == mq
                    % Overlap occurs, subtract volume of smaller member
                    xk1 = NC(SortedCA(k,1),1); xk2 = NC(SortedCA(k,2),1);
                    yk1 = NC(SortedCA(k,1),2); yk2 = NC(SortedCA(k,2),2);
                    Lk = sqrt(((xk2-xk1)^2)+((yk2-yk1)^2));
                    xq1 = NC(SortedCA(q,1),1); xq2 = NC(SortedCA(q,2),1);
                    yq1 = NC(SortedCA(q,1),2); yq2 = NC(SortedCA(q,2),2);
                    Lq = sqrt(((xq2-xq1)^2)+((yq2-yq1)^2));
                    if Lk >= Lq
                        tTV = tTV-(0.5*pi*(r^2)*Lq);
                        disp(tTV);
                    elseif Lk < Lq
                        tTV = tTV-(0.5*pi*(r^2)*Lk);
                        disp(tTV);
                    end
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
                   % Overlap occurs, subtract volume of smaller member
                    xk1 = NC(SortedCA(k,1),1); xk2 = NC(SortedCA(k,2),1);
                    yk1 = NC(SortedCA(k,1),2); yk2 = NC(SortedCA(k,2),2);
                    Lk = sqrt(((xk2-xk1)^2)+((yk2-yk1)^2));
                    xq1 = NC(SortedCA(q,1),1); xq2 = NC(SortedCA(q,2),1);
                    yq1 = NC(SortedCA(q,1),2); yq2 = NC(SortedCA(q,2),2);
                    Lq = sqrt(((xq2-xq1)^2)+((yq2-yq1)^2));
                    if Lk >= Lq
                        tTV = tTV-(0.5*pi*(r^2)*Lq);
                        disp(tTV);
                    elseif Lk < Lq
                        tTV = tTV-(0.5*pi*(r^2)*Lk);
                        disp(tTV);
                    end
                end
            end
        end
    end
    %}
end

% FUNCTION TO DETERMINE PRESENCE OF INTERSECTION (FOR CONSTRAINT #1)
function [intersect,alpha] = findLineSegIntersection(p1,q1,p2,q2)
    if (findOrientation(p1,q1,p2) ~= findOrientation(p1,q1,q2))&&...
            (findOrientation(p2,q2,p1) ~= findOrientation(p2,q2,q1)) 
        if isequal(p1,p2) || isequal(q1,q2) || ...
                isequal(p1,q2) || isequal(q1,p2)
            intersect = false; alpha = 0;
        else
            intersect = true;
            v1 = [q1(1),q1(2),0] - [p1(1),p1(2),0];
            v2 = [q2(1),q2(2),0] - [p2(1),p2(2),0];
            alpha = acos(max(min(dot(v1,v2)/(norm(v1)*norm(v2)),1),-1));
        end
    else
        intersect = false; alpha = 0;
    end
end

% FUNCTION TO CALCULATE ORIENTATION FROM 3 POINTS (FOR CONSTRAINT #1)
function orientation = findOrientation(p,q,r)
    val = ((q(2)-p(2))*(r(1)-q(1)))-((q(1)-p(1))*(r(2)-q(2)));
    if val == 0
        orientation = 0;
    elseif val > 0
        orientation = 1;
    else
        orientation = 2;
    end
end

% FUNCTION TO SUBTRACT VOLUME OVERLAP AT NODES
function tTV = subNodOLVol(NC,CA,tTV,rvar)
    for z = 1:1:size(NC,1)
        % Isolate members originating/ending at the current node
        indione = CA(:,1) == z; inditwo = CA(:,2) == z;
        mCAone = CA(indione,:); mCAtwo = CA(inditwo,:);
        mCA = ...
            [setdiff(mCAone,mCAtwo,'rows');setdiff(mCAtwo,mCAone,'rows')];
        
        % Loop through each pair of connecting members at the current node
        for i = 1:1:size(mCA,1)
            for j = 1:1:size(mCA,1)
                if i ~= j
                    % Find the sweep angle between the connecting members
                    p1 = [NC(mCA(i,1),1),NC(mCA(i,1),2)];
                    q1 = [NC(mCA(i,2),1),NC(mCA(i,2),2)];
                    p2 = [NC(mCA(j,1),1),NC(mCA(j,1),2)];
                    q2 = [NC(mCA(j,2),1),NC(mCA(j,2),2)];
                    v1 = [q1(1),q1(2),0] - [p1(1),p1(2),0];
                    v2 = [q2(1),q2(2),0] - [p2(1),p2(2),0];
                    alpha = ...
                       acos(max(min(dot(v1,v2)/(norm(v1)*norm(v2)),1),-1));
                    theta = pi - alpha;
                    
                    % Check if members overlap
                    if (abs(theta-pi) < (10^-5)) || (abs(theta) < (10^-5))
                        continue
                    else
                        % Find the volume of each overlap sphere sector,
                        % subtract from the total truss volume
                        frac = theta/(2*pi);
                        [~,idx1] = ismember(mCA(i,:),CA,'rows');
                        [~,idx2] = ismember(mCA(j,:),CA,'rows');
                        r1 = rvar(idx1); r2 = rvar(idx2);
                        avgrad = mean([r1,r2]);
                        VOL = ((4*pi)/3)*(avgrad^3)*frac;
                        tTV = tTV - VOL;
                        %disp(tTV);
                    end
                end
            end
        end
    end
end