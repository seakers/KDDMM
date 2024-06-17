% -----------------------------------------------------------------
% ANSYS APDL 2D NxN BEAM CODE, with Periodic Boundary Conditions
% This code considers only NxN 2D designs, with constant member radii
% This code calculates metamaterial properties for a single unit cell
% All lengths are in [m], all stresses and moduli are in [Pa] 
% Each node has two degrees of freedom (DOF), x and y.  The (node number*2) 
%    gives the number for that node's yDOF, and (yDOF - 1) is the xDOF
% -----------------------------------------------------------------
% This code uses MATLAB to call ANSYS APDL for FEM analysis.
% First, MATLAB writes the design parameters into input file 
% (e.g. para_in.txt)
% Then, MATLAB calls ANSYS APDL to excute the APDL file 
% The APDL file reads the parameters from the input file and writes the 
% analysis result to the output file (e.g. para_out.txt)
% Finally, MATLAB reads the results from the output file
% -----------------------------------------------------------------
% CREDITS FOR ORIGINAL CODE PACKAGE:
% https://github.com/zhandawei/Call_ANSYS_in_MATLAB
% Dawei Zhan
% zhandawei{at}swjtu.edu.cn
% -----------------------------------------------------------------
% Modified/adapted for KD3M2 project by Srikar Srivatsa
% -----------------------------------------------------------------
% Sample values for Case 1 (for testing script):
%{
clc;   
close all; 
clear;
sidenum = 3;
sel = 0.01; 
E = 1816200;
r = 250*(10^-6);
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;2,5;3,5;4,5;5,6;5,7;5,8;5,9];
%}

function C = Beam_2D_NxN_PBC(sel,sidenum,r,E,CA)
    % Initialize parameters
    C = [0,0,0;0,0,0;0,0,0];
    
    % Generate nodal grid
    NC = generateNC(sel,sidenum);
    
    % Determine edge members
    [CA,edgelog] = edgeVals(CA,sidenum);
    
    % Calculate cross-sectional properties for edge members
    [csPropsA,csPropsB] = calcCSProps(r);
    
    % Determine nodal pairs for each constraint equation
    [NP,numLRCE,numTBCE] = findNodalPairs(sidenum,CA,NC);
    
    % Iterate through once for each strain component
    for y = 1:1:3
        % Define strain vector: [e11, e22, e12]'
        strainvec = [0;0;0];

        % Set that component equal to a dummy value (0.01 strain), set all 
        % other values to zero
        strainvec(y) = 0.001; 
        strainvec(3) = strainvec(3)*2;
        
        % Populate Input Vector
        qsel = sel./(10^-3); qr = r./(10^-6); 
        qcsPropsA = csPropsA./(10^-6); qcsPropsB = csPropsB./(10^-15);
        input_vec=[qsel;qr;E;size(CA,1);y;sidenum;numLRCE;numTBCE;...
          qcsPropsA';qcsPropsB';edgelog;NP(:,1);NP(:,2);CA(:,1);CA(:,2)];

        % Define paths and input/output files 
        ANSYS_path=...
 'C:\Program Files\ANSYS Inc\v211\ansys\bin\winx64\ANSYS211';
        APDL_name='C:\Users\Public\APDL_Runs\Truss2D_NxN.txt';
        input_file_name='C:\Users\Public\APDL_Runs\para_in.txt';
        output_file_name='C:\Users\Public\APDL_Runs\para_out.txt';
        ANSYS_path=strcat('"',ANSYS_path,'"');
        disp('Current stage: ');disp(y);

        % Write input vector to para_in.txt
        %writematrix(input_vec,input_file_name);
        dlmwrite(input_file_name,input_vec,'delimiter',' ',...
           'precision','%8.5f','newline','pc');
        % Call ANSYS APDL, perform FEA
        status = system(sprintf('%s -b -p aa_r -i %s -o out.txt -np 4',...
            ANSYS_path,APDL_name));
        disp('Simulation status: ');disp(status);
        % Read the results from para_out.txt
        output_vec=load(output_file_name)';
        F_x = output_vec(1);
        F_y = output_vec(2);
        F_xy = output_vec(3);
        FBank(y,:) = [F_x,F_y,F_xy];

        % Calculate stress vector
        stressvec = (1/(sel*2*r)).*[F_x;F_y;F_xy];

        % Use strain and stress vectors to solve for the corresponding row 
        % of the C matrix
        Cdummy = stressvec/strainvec;
        C(:,y) = Cdummy(:,y);
    end
end

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

% FUNCTION TO CALCULATE EDGE-RELATED INPUTS
function [CA,edgelog] = edgeVals(CA,sidenum)
       
    % Identify four sets of edge nodes
    edge1nodes = 1:1:sidenum;
    edge2nodes = ((sidenum^2)-sidenum+1):1:(sidenum^2);
    edge3nodes = 1:sidenum:((sidenum^2)-(sidenum)+1);
    edge4nodes = sidenum:sidenum:(sidenum^2);
             
    % Identify edge members and edge types
    edge1connA = ismember(CA(:,1),edge1nodes);
    edge1connB = ismember(CA(:,2),edge1nodes);
    edge1log = (edge1connA & edge1connB); 
    edge1type = (edge1connA & edge1connB); 
    edge2connA = ismember(CA(:,1),edge2nodes);
    edge2connB = ismember(CA(:,2),edge2nodes);
    edge2log = (edge2connA & edge2connB); 
    edge2type = 2.*(edge2connA & edge2connB);
    edge3connA = ismember(CA(:,1),edge3nodes);
    edge3connB = ismember(CA(:,2),edge3nodes);
    edge3log = (edge3connA & edge3connB); 
    edge3type = 3.*(edge3connA & edge3connB);
    edge4connA = ismember(CA(:,1),edge4nodes);
    edge4connB = ismember(CA(:,2),edge4nodes);
    edge4log = (edge4connA & edge4connB); 
    edge4type = 4.*(edge4connA & edge4connB);
    
    edgelog = edge1log + edge2log + edge3log + edge4log;
    edgetype = edge1type + edge2type + edge3type + edge4type;
    
    % Modify order of edge members depending on type
    for i = 1:1:size(CA,1)
        if (edgetype(i) == 1) || (edgetype(i) == 4)
            if CA(i,2) > CA(i,1)
                newrow = [CA(i,2),CA(i,1)];
                CA(i,:) = newrow;
            end
        elseif (edgetype(i) == 2) || (edgetype(i) == 3)
            if CA(i,1) > CA(i,2)
                newrow = [CA(i,2),CA(i,1)];
                CA(i,:) = newrow;
            end
        end
    end
end

% FUNCTION TO CALCULATE CROSS-SECTIONAL PROPERTIES OF EDGE MEMBERS
function [csPropsA,csPropsB] = calcCSProps(r)
    A = 0.5*pi*(r^2); % Area of section
    Iyy = (pi*(r^4))/8; % Moment of inertia about the y axis
    Izz = (pi*(r^4))/8; % Moment of inertia about the z axis
    Iw = 0; % Warping constant
    J = (pi*(r^4))/4; % Torsional constant
    CGy = 0; % y coordinate of centroid
    CGz = (4*r)/(3*pi); % z coordinate of centroid
    Iyz = A*CGy*CGz; % Product of inertia
    SHy = 0; % y coordinate of shear center
    SHz = (4*r)/(3*pi); % z coordinate of shear center
    TKz = r; % Thickness along Z axis (maximum height)
    TKy = r/2; % Thickness along Y axis (maximum width)
    csPropsA = [A,Iyz,Iw,CGy,CGz,SHy,SHz,TKz,TKy];
    csPropsB = [Iyy,Izz,J];
end

% FUNCTION TO DETERMINE NODAL PAIRS FOR CONSTRAINTS EQUATIONS
function [NP,numLRCE,numTBCE] = findNodalPairs(sidenum,CA,NC)
    % Add up counters based on nodal connectivities 
    edgein = 0.5:1:(0.5+size(NC,1));
    [N,~] = histcounts(CA,edgein);
    
    % Determine vertical-edge nodal pairs
    NP = []; numLRCE = 0;
    for i = 1:1:sidenum
        if (N(i) == 0) && (N(i+(sidenum*(sidenum-1))) == 0)
            % Both nodes are unoccupied, no constraint equation required
        elseif (N(i) == 0) && (N(i+(sidenum*(sidenum-1))) ~= 0)
            % Left-edge node is unoccupied
            index = 0; counter = 0;
            while index == 0
                lower = 1 + counter; 
                upper = sidenum + counter;
                Nlocal = N(lower:upper);
                index = find(Nlocal,1,'first');
                counter = counter + sidenum;
            end
            NP = [NP;(index+(counter-sidenum)),(i+(sidenum*(sidenum-1)))];
            numLRCE = numLRCE + 1;
        elseif (N(i) ~= 0) && (N(i+(sidenum*(sidenum-1))) == 0)
            % Right-edge node is unoccupied
            index = 0; counter = 0;
            while index == 0
                lower = ((sidenum^2)-sidenum+1) - counter; 
                upper = (sidenum^2) - counter;
                Nlocal = N(lower:upper);
                index = find(Nlocal,1,'first');
                counter = counter + sidenum;
            end
            NP = [NP;i,(index+((sidenum^2)-sidenum)-(counter-sidenum))];
            numLRCE = numLRCE + 1;
        else
            % Both nodes are occupied, substitute nodes not required
            NP = [NP;i,(i+(sidenum*(sidenum-1)))];
            numLRCE = numLRCE + 1;
        end
    end
    
    % Determine horizontal-edge nodal pairs
    numTBCE = 0;
    for i = 1:sidenum:((sidenum^2)-sidenum+1)
        if (N(i) == 0) && (N(i+sidenum-1) == 0)
            % Both nodes are unoccupied, no constraint equation required
        elseif (N(i) == 0) && (N(i+sidenum-1) ~= 0)
            % Bottom-edge node is unoccupied
            index = 0; counter = 0;
            while index == 0
                lower = 1 + counter; 
                upper = (sidenum^2)-sidenum+1 + counter;
                Nlocal = N(lower:upper);
                index = find(Nlocal,1,'first');
                counter = counter + 1;
            end
            NP = [NP;((((index-1)*sidenum)+1)+(counter-1)),(i+sidenum-1)];
            numTBCE = numTBCE + 1;
        elseif (N(i) ~= 0) && (N(i+sidenum-1) == 0)
            % Top-edge node is unoccupied
            index = 0; counter = 0;
            while index == 0
                lower = sidenum - counter; 
                upper = (sidenum^2) - counter;
                Nlocal = N(lower:upper);
                index = find(Nlocal,1,'first');
                counter = counter + 1;
            end
            NP = [NP;i,((index*sidenum)-(counter-1))];
            numTBCE = numTBCE + 1;
        else
            % Both nodes are occupied, substitute nodes not required
            NP = [NP;i,(i+sidenum-1)];
            numTBCE = numTBCE + 1;
        end
    end
end



