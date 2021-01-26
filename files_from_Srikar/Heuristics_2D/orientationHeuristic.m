% 2D Orientation Heuristic
% This function evaluates the orientation of members in a 2D NxN design,
% and assigns a score based on whether the average orientation of members
% is biased towards either the vertical or horizontal directions.  The
% target angle is based off of a logistic function
% -----------------------------------------------------------------
% Sample values for Case 1 (for testing script):
%{
target = 1;
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;2,5;3,5;4,5;5,6;5,7;5,8;5,9];
NC = [0,0;0,0.025;0,0.05;0.025,0;0.025,0.025;
      0.025,0.05;0.05,0;0.05,0.025;0.05,0.05];
%}

function [orientationScore,avgAngle] = orientationHeuristic(NC,CA,target)
    % Finding angles of members relative to horizontal
    x1 = NC(CA(:,1),1); x2 = NC(CA(:,2),1);
    y1 = NC(CA(:,1),2); y2 = NC(CA(:,2),2);
    L = sqrt(((x2-x1).^2)+((y2-y1).^2));
    angles = abs(acos((x2-x1)./L));
    
    % Find average angle
    avgAngle = rad2deg(mean(angles));
    
    % Find target angle
    if target < 1
        ntarget = 1/target;
    else
        ntarget = target;
    end
    k = 5;
    tA = 90./(1+exp(-k.*(ntarget-1)));
    if target < 1
        targetAngle = 90 - tA;
    else
        targetAngle = tA;
    end
    
    % Find deviation of average angle from target angle
    deviation = avgAngle - targetAngle;
    
    % Score deviation of average angle
    orientationScore = 1 - abs(0.01.*deviation);
    if orientationScore < 0
        orientationScore = 0;
    end
end