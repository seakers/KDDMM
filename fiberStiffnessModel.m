% Fiber Stiffness-based model ("Simple Model")
% This model is intended to find C-matrix (material stiffness) values based
%    approximations of individual truss stiffnesses from treating them as
%    fibers at various orientations
% All lengths are in [m], all stresses and moduli are in [Pa] 
% The model uses any 2D NxN nodal grid
% Each node has two degrees of freedom (DOF), x and y.  The (node number*2) 
%    gives the number for that node's yDOF, and (yDOF - 1) is the xDOF
function [C11,C22] = fiberStiffnessModel(sel,r,E,CA,sidenum,nucFac)
    %disp(strcat('sel = ',num2str(sel)))
    %disp(strcat('r = ',num2str(r)))
    %disp(strcat('E = ',num2str(E)))
    %disp(strcat('sidenum = ',num2str(sidenum)))
    %disp(strcat('nucFac = ',num2str(nucFac)))
    %disp(CA)
    % Generate nodal grid
    NC = generateNC(sel,sidenum);
    %disp(NC)
    
    % Find volume fraction
    volFrac = calcVF(NC,CA,r,sel);
    %disp(volFrac)
    
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
    %disp(K)
    
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

%%
% sel = 10e-3;
% r = 250e-6;
% E = 1.8162e6;
% CA = [1,2;2,3;1,4;1,5;2,5;3,5;3,6;4,5;5,6;4,7;5,7;5,8;5,9;6,9;7,8;8,9];
% sidenum = 3;
% nucFac = 3;
% [C11,C22] = fiberStiffnessModel(sel,r,E,CA,sidenum,nucFac)