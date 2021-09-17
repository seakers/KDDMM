% Fiber Stiffness-based model ("Simple Model")
% This model is intended to find C-matrix (material stiffness) values based
%    approximations of individual truss stiffnesses from treating them as
%    fibers at various orientations
% All lengths are in [m], all stresses and moduli are in [Pa] 
% The model uses any 2D NxN nodal grid
% Each node has two degrees of freedom (DOF), x and y.  The (node number*2) 
%    gives the number for that node's yDOF, and (yDOF - 1) is the xDOF
% -----------------------------------------------------------------
% Sample values to test code (copy-paste into command window):
%{
clc;    
close all; 
clear;
nucFac = 1; 
sel = 0.05; 
E = 10000;
% Case 1, 1 unit cell (3x3 grid)
CA = [1,2;2,3;1,4;1,5;2,5;3,5;3,6;4,5;5,6;4,7;5,7;5,8;5,9;6,9;7,8;8,9];
rvar = (50*(10^-6)).*ones(1,size(CA,1)); rvar(1) = (55*(10^-6));
%}
function [C11,C22,volFrac] = ...
                    fiberStiffnessModel_rVar_V2(sel,rvar,E,CA,nucFac)
    % Calculated Inputs
    sidenum = (2*nucFac) + 1;
    Avar = pi.*(rvar.^2); % Cross-sectional areas of truss members

    % Generate nodal grid
    NC = generateNC(sel,sidenum);
    
    % Find volume fraction
    volFrac = calcVF(NC,CA,rvar,sel,sidenum);
    
    % Calculating C-matrix values
    C11 = fiberCalc(NC,CA,volFrac,Avar,E,1)/nucFac;
    C22 = fiberCalc(NC,CA,volFrac,Avar,E,2)/nucFac;
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

% FUNCTION TO CALCULATE C-MATRIX VALUES VIA FIBER METHOD
function Cval = fiberCalc(NC,CA,volFrac,Avar,E,dir)    
    % Find effective structural stiffness
    K = E*volFrac;
    
    % Find volume-corrected sum of cosines for all fibers 
    cVsum = 0;
    Vsum = 0;
    for i = 1:size(CA,1)
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        if dir == 1
            c=(x2-x1)/L; 
        elseif dir == 2
            c=(y2-y1)/L; 
        end
        cVsum = cVsum + (L*(Avar(i))*(c^4));
        Vsum = Vsum + (L*Avar(i));
    end
    
    % Find desired C-value
    Cval = (K*cVsum)/Vsum;
end

% FUNCTION TO CALCULATE VOLUME FRACTION 
function volFrac = calcVF(NC,CA,rvar,sel,sidenum)
    totalTrussVol = 0;
    for i = 1:size(CA,1)
        % Finding element length from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        % Adding current element to total volume of trusses
        totalTrussVol = totalTrussVol + (L*pi*(rvar(i)^2));
    end
    
    % Finding average side "thickness" due to differing element radii
    horizrads = [];
    for i = 1:1:size(CA,1)
        if ((CA(i,1) + sidenum) == CA(i,2)) && (NC(CA(i,1),2) == sel)
            horizrads = [horizrads,rvar(i)];
        end
    end
    vertrads = [];
    for i = 1:1:size(CA,1)
        if ((CA(i,1) + 1) == CA(i,2)) && (NC(CA(i,1),1) == sel)
            vertrads = [vertrads,rvar(i)];
        end
    end
    thick = mean([mean(horizrads),mean(vertrads)]);
    
    % Calculating volume fraction (using a solid square with 2*(avg 
    %   thickness) as a baseline)
    volFrac = totalTrussVol/(2*thick*(sel^2)); 
end