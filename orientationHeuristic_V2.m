% 2D Orientation Heuristic
% This function evaluates the orientation of members in a 2D NxN design,
% and assigns a score based on whether the average orientation of members
% is biased towards either the vertical or horizontal directions.  The
% target angle is based off of an arctan function
% -----------------------------------------------------------------
% Sample values for Case 1 (for testing script):
%{
target = 1;
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;2,5;3,5;4,5;5,6;5,7;5,8;5,9];
NC = [0,0;0,0.025;0,0.05;0.025,0;0.025,0.025;
      0.025,0.05;0.05,0;0.05,0.025;0.05,0.05];
%}

function [orientationScore,avgAngle] = ...
                     orientationHeuristic_V2(NC,CA,target)
    % Finding angles of members relative to horizontal
    SortedCA = sortrows((sort(CA'))');
    x1 = NC(CA(:,1),1); x2 = NC(CA(:,2),1);
    y1 = NC(CA(:,1),2); y2 = NC(CA(:,2),2);
    angles = atan((y2-y1)./(x2-x1));
    for i = 1:1:length(angles)
        if angles(i)<0
            angles(i) = abs(angles(i));
        end
    end
    
    % Find average angle
    avgAngle = rad2deg(mean(angles));
    
    % Find target angle
    targetAngle = rad2deg(atan(target));
    
    % Find deviation of average angle from target angle
    deviation = avgAngle - targetAngle;
    
    % Score deviation of average angle
    orientationScore = 1 - abs(0.01.*round(deviation,-1));
    if orientationScore < 0
        orientationScore = 0;
    end
end
