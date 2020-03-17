% Checking Script for Constraint Checker V2 (based on 2D NxN Truss Code)

% Initialization
clc;    
close all; 
clear;

% Input Constants
sidenum = 3; % Number of nodes on each size of square grid
nucFac = 2; % factor to indicate how many unit cells are in this model
            %   (so nucFac = 1 indicates the entire 5x5 is one unit cell,
            %    and nucFac = 2 indicates that the 5x5 contains four
            %    3x3 unit cells inside of it
sel = 0.05; % unit cube side length of 5 cm, truss length of 2.5 cm
r = 50*(10^-6); % Radius of 50 micrometers for cross-sectional 
                %   area of (assumed circular) truss members
A = pi*(r^2); % Cross-sectional area of truss member
E = 10000; % Young's Modulus of 10000 Pa (polymeric material)

% Input trial connectivity arrays
% Case 1, 1 unit cell
CAone = ...
    [1,2;2,3;1,4;1,5;2,5;3,5;3,6;4,5;5,6;4,7;5,7;5,8;5,9;6,9;7,8;8,9];
% Case 2, 1 unit cell
CAtwo = ...
    [1,2;2,3;1,4;2,4;2,5;2,6;3,6;4,5;5,6;4,7;4,8;5,8;6,8;6,9;7,8;8,9];
% Combination of case 1 and 2 (crossing diagonal members)
CAthree = ...
    [1,2;2,3;1,4;2,4;2,5;2,6;3,6;4,5;5,6;4,7;4,8;5,8;6,8;6,9;7,8;8,9;...
     1,5;3,5;5,7;5,9];

% Generate vector with nodal coordinates
NC = generateNC(sel,sidenum); 

% Check constraints for all trial truss designs
ver_one = constChecker_V2(NC,CAone);
ver_two = constChecker_V2(NC,CAtwo);
ver_three = constChecker_V2(NC,CAthree);

%----------%
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

% FUNCTION FOR CONSTRAINT CHECKER
% This function is intended to check that a given truss design is
% legitimate (not necessarily stable, but legitimate within the standard
% nodal framework and related assumptions)
function constVerified = constChecker_V2(NC,CA)
    constVerified = true;

    % FIRST CONSTRAINT: members only intersect at nodes (no crossing)
    % (source: https://www.geeksforgeeks.org/given-a-set-of-line-segments
    % -find-if-any-two-segments-intersect/)
    
    % Sort points from left to right by x-position
    SortedCA = sortrows(CA);
    
    % Develop 4xM matrix of line segment endpoint coordinates, where M is 
    %   the number of truss members.  Each row of format (x1,y1,x2,y2),
    %   where point 1 is leftmost, point 2 is rightmost
    PosA = [NC(SortedCA(:,1),1),NC(SortedCA(:,1),2),...
            NC(SortedCA(:,2),1),NC(SortedCA(:,2),2)];
        
    for i = 1:1:size(PosA,1)
        for j = 1:1:size(PosA,1)
            intersect = findLineSegIntersection([PosA(i,1),PosA(i,2)],...
                        [PosA(i,3),PosA(i,4)],[PosA(j,1),PosA(j,2)],...
                        [PosA(j,3),PosA(j,4)]);
            if intersect == true
                constVerified = false;
                D = ['The element from (',num2str(PosA(i,1)),',',...
                     num2str(PosA(i,2)),') to (',num2str(PosA(i,3)),...
                     ',',num2str(PosA(i,4)),') intersects with the'...
                     ' element from (',num2str(PosA(j,1)),',',...
                     num2str(PosA(j,2)),') to (',num2str(PosA(j,3)),',',...
                     num2str(PosA(j,4)),')'];
                disp(D);
                return
            end
        end
    end
end

% FUNCTION TO DETERMINE PRESENCE OF INTERSECTION
% (source: https://www.geeksforgeeks.org/check-if-two-given-line-segments
% -intersect/)
% This boolean function determines whether two line segments intersect,
% given their endpoints as inputs
function intersect = findLineSegIntersection(p1,q1,p2,q2)
    if (findOrientation(p1,q1,p2) ~= findOrientation(p1,q1,q2))&&...
            (findOrientation(p2,q2,p1) ~= findOrientation(p2,q2,q1)) 
        if isequal(p1,p2) || isequal(q1,q2) || ...
                isequal(p1,q2) || isequal(q1,p2)
            intersect = false; 
        else
            intersect = true;
        end
    else
        intersect = false;
    end
end

% FUNCTION TO CALCULATE ORIENTATION FROM 3 POINTS
% (source: https://www.geeksforgeeks.org/orientation-3-ordered-points/)
% This function finds the orientation of an ordered triplet (p, q, r)
% The function returns one of three following values:
% 0 --> p, q and r are colinear 
% 1 --> Clockwise 
% 2 --> Counterclockwise 
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