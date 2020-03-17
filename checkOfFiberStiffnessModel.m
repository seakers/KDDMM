% Checking Script for Fiber Stiffness Model (based on 2D NxN Truss Code)

% Initialization
clc;    
close all; 
clear;

% Input Constants (3x3)
sidenum_3x3 = 3; % Number of nodes on each size of square grid
nucFac_3x3 = 1; % factor to indicate how many unit cells are in this model
            %   (so nucFac = 1 indicates the entire 5x5 is one unit cell,
            %    and nucFac = 2 indicates that the 5x5 contains four
            %    3x3 unit cells inside of it
sel = 0.05; % unit cube side length of 5 cm, truss length of 2.5 cm
r = 50*(10^-6); % Radius of 50 micrometers for cross-sectional 
                %   area of (assumed circular) truss members
E = 10000; % Young's Modulus of 10000 Pa (polymeric material)

% Input trial connectivity arrays (3x3)
% Case 1, 1 unit cell
CAone = ...
    [1,2;2,3;1,4;1,5;2,5;3,5;3,6;4,5;5,6;4,7;5,7;5,8;5,9;6,9;7,8;8,9];
% Case 2, 1 unit cell
CAtwo = ...
    [1,2;2,3;1,4;2,4;2,5;2,6;3,6;4,5;5,6;4,7;4,8;5,8;6,8;6,9;7,8;8,9];
 
% Check constraints for all trial truss designs (3x3)
ver_one = fiberStiffnessModel(sel,r,E,CAone,sidenum_3x3,nucFac_3x3);
ver_two = fiberStiffnessModel(sel,r,E,CAtwo,sidenum_3x3,nucFac_3x3);

% Input Constants (5x5)
sidenum_5x5 = 5; % Number of nodes on each size of square grid
nucFac_5x5 = 2; % factor to indicate how many unit cells are in this model
            
% Input trial connectivity arrays (5x5)
% Case 1, 4 unit cells
CAthree = [1,2;2,3;3,4;4,5;5,10;10,15;15,20;20,25;
           25,24;24,23;23,22;22,21;16,21;11,16;6,11;1,6;
           3,8;8,13;13,18;18,23;11,12;12,13;13,14;14,15;
           1,7;2,7;3,7;6,7;7,8;7,11;7,12;7,13;
           3,9;4,9;5,9;8,9;9,10;9,13;9,14;9,15;
           11,17;12,17;13,17;16,17;17,18;17,21;17,22;17,23;
           13,19;14,19;15,19;18,19;19,20;19,23;19,24;19,25];
% Case 2, 4 unit cells
CAfour = [1,2;2,3;3,4;4,5;5,10;10,15;15,20;20,25;
          25,24;24,23;23,22;22,21;16,21;11,16;6,11;1,6;
          3,8;8,13;13,18;18,23;11,12;12,13;13,14;14,15;
          2,6;2,7;2,8;6,7;7,8;6,12;7,12;8,12;
          4,8;4,9;4,10;8,9;9,10;8,14;9,14;10,14;
          12,16;12,17;12,18;16,17;17,18;16,22;17,22;18,22;
          14,18;14,19;14,20;18,19;19,20;18,24;19,24;20,24];
 
% Check constraints for all trial truss designs (5x5)
ver_three = fiberStiffnessModel(sel,r,E,CAthree,sidenum_5x5,nucFac_5x5);
ver_four = fiberStiffnessModel(sel,r,E,CAfour,sidenum_5x5,nucFac_5x5);

%----------------------%
% FUNCTION FOR FIBER STIFFNESS MODEL
function [C11,C22] = fiberStiffnessModel(sel,r,E,CA,sidenum,nucFac)
    % Generate nodal grid
    NC = generateNC(sel,sidenum);
    
    % Find volume fraction
    volFrac = calcVF(NC,CA,r,sel);
    
    % Calculating C-matrix values
    C11 = fiberCalc(volFrac,NC,CA,E,1)/nucFac;
    C22 = fiberCalc(volFrac,NC,CA,E,2)/nucFac;

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

% FUNCTION TO CALCULATE VOLUME FRACTION
function volFrac = calcVF(NC,CA,r,sel)
    totalTrussVol = 0;
    for i = 1:size(CA,1)
        % Finding element length from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        % Adding current element to total volume of trusses
        totalTrussVol = totalTrussVol + (L*pi*(r^2));
    end
    % Calculating volume fraction (using a solid square with 2*r thickness
    %   as a baseline)
    volFrac = totalTrussVol/(2*r*(sel^2));
end

% FUNCTION TO CALCULATE C-MATRIX VALUES VIA FIBER METHOD
function Cval = fiberCalc(volFrac,NC,CA,E,dir)
    % Find effective structural stiffness
    K = E*volFrac;
    
    % Find length-corrected sum of cosines for all fibers 
    cLsum = 0;
    Lsum = 0;
    for i = 1:size(CA,1)
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        if dir == 1
            c=(x2-x1)/L; 
        elseif dir == 2
            c=(y2-y1)/L; 
        end
        cLsum = cLsum + (L*(c^4));
        Lsum = Lsum + L;
    end
    
    % Find desired C-value
    Cval = (K*cLsum)/Lsum;
    
end