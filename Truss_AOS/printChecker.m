% 2D Printability Checker
% This function evaluates whether a 2D NxN design can be realistically
% printed, retaining its structural and geometric integrity as intended
% -----------------------------------------------------------------
% Sample values for Case 1 (for testing script):
%{
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;2,5;3,5;4,5;5,6;5,7;5,8;5,9];
NC = [0,0;0,0.025;0,0.05;0.025,0;0.025,0.025;
      0.025,0.05;0.05,0;0.05,0.025;0.05,0.05];
rvar = (50*(10^-6)).*ones(1,size(CA,1));
%}
function [printScore,printBool] = printChecker(NC,CA,rvar)
    % Initialize printScore, printBool
    printBool = true;
    printScore = 1;

    % Length Allowance
    x1 = NC(CA(:,1),1); x2 = NC(CA(:,2),1);
    y1 = NC(CA(:,1),2); y2 = NC(CA(:,2),2);
    L = sqrt(((x2-x1).^2)+((y2-y1).^2));
    for i = 1:1:size(L,1)
        % Bridge Criterion
        if L(i) > 0.006
            printScore = printScore - 0.1;
            if printScore < 0.1
                printBool = false;
                return
            end
        end
        % Aspect Ratio Criterion
        if (4*rvar(i)) > L(i)
            printScore = printScore - 0.1;
            if printScore < 0.1
                printBool = false;
                return
            end
        end
    end
    
    % Angle Restrictions
    % Loop through each node
    for i = 1:1:size(NC,1)
        % Isolate elements originating/ending at a given node
        indione = CA(:,1) == i; inditwo = CA(:,2) == i;
        mCAone = CA(indione,:); mCAtwo = CA(inditwo,:);
        mCA = ...
            [setdiff(mCAone,mCAtwo,'rows');setdiff(mCAtwo,mCAone,'rows')];
        for q = 1:1:size(mCA,1)
            if mCA(q,2) == i
                mCA(q,:) = [mCA(q,2),mCA(q,1)];
            end
        end
        
        % Reorder mCA by elements' angles with the horizontal so that
        % subsequent elements in mCA are physically adjacent
        mx1 = NC(mCA(:,1),1); mx2 = NC(mCA(:,2),1);
        my1 = NC(mCA(:,1),2); my2 = NC(mCA(:,2),2);
        mL = sqrt(((mx2-mx1).^2)+((my2-my1).^2));
        phis = rad2deg(acos((mx2-mx1)./mL)); phis2 = phis;
        for t = 1:1:size(phis,1)
            if (mx2(t)-mx1(t)) >= 0 && (my2(t)-my1(t)) >= 0
                % Angle in 1st quadrant; do nothing
            elseif (mx2(t)-mx1(t)) < 0 && (my2(t)-my1(t)) >= 0
                % Angle in second quadrant; do nothing
            elseif (mx2(t)-mx1(t)) < 0 && (my2(t)-my1(t)) < 0
                % Angle in third quadrant; add 90 degrees
                phis(t) = phis(t) + 90;
            elseif (mx2(t)-mx1(t)) == 0 && (my2(t)-my1(t)) < 0
                % Angle is 270 degrees; add 180 degrees
                phis(t) = phis(t) + 180;
            elseif (mx2(t)-mx1(t)) > 0 && (my2(t)-my1(t)) < 0
                % Angle in fourth quadrant; add 270 degrees
                phis(t) = phis(t) + 270;
            end
        end
        [~,indices] = sort(phis);
        msCA = mCA(indices,:);
        
        % For each pair of adjacent elements at the current node (in mCA),
        % find and evaluate the angle inbetween
        for j = 1:1:size(msCA,1)
            % Find slopes of adjacent line segments
            m1 = (NC(msCA(j,2),2)-NC(msCA(j,1),2))/...
                 (NC(msCA(j,2),1)-NC(msCA(j,1),1));
            if j == size(msCA,1)
                m2 = (NC(msCA(1,2),2)-NC(msCA(1,1),2))/...
                     (NC(msCA(1,2),1)-NC(msCA(1,1),1));
            else
                m2 = (NC(msCA(j+1,2),2)-NC(msCA(j+1,1),2))/...
                     (NC(msCA(j+1,2),1)-NC(msCA(j+1,1),1));      
            end
            
            % Find enclosed angle
            theta = atan(abs((m1-m2)/(1+(m1*m2))));
            
            % Judge angle
            if rad2deg(theta) < 40
                printScore = printScore - 0.1;
                if printScore < 0.1
                    printBool = false;
                    return
                end
            end
        end
    end
    
    % Assign printBool
    if printScore < 1
        printBool = false;
    end
end
