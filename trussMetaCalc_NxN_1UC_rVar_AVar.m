% -----------------------------------------------------------------
% Function for calculating metamaterial properties for a 2D NxN truss with
% individually-adjustable element radii, for a single unit cell
% -----------------------------------------------------------------
% All lengths are in [m], all stresses and moduli are in [Pa] 
% Each node has two degrees of freedom (DOF), x and y.  The (node number*2) 
%    gives the number for that node's yDOF, and (yDOF - 1) is the xDOF
% -----------------------------------------------------------------
% This function accounts for repetition of unit cells by considering only
% half the cross-sectional area of edge members
% -----------------------------------------------------------------
% Sample values to test code (copy-paste into command window):
%{
clc;    
close all; 
clear;
nucFac = 3; 
sel = 0.01; 
E = 1816200;
% Case 1, 1 unit cell (3x3 grid)
CA = [1,2;2,3;1,4;1,5;2,5;3,5;3,6;4,5;5,6;4,7;5,7;5,8;5,9;6,9;7,8;8,9];
rvar = (250*(10^-6)).*ones(1,size(CA,1));
%}
function [C,volFrac] = ...
             trussMetaCalc_NxN_1UC_rVar_AVar(sidenum,sel,rvar,E,CA)
    
    % Generate vector with nodal coordinates
    NC = generateNC(sel,sidenum);
    
    % Calculate Avar & modify for edge members
    Avar = pi.*(rvar.^2); % Cross-sectional areas of truss members
    Avar = modifyAreas(Avar,CA,NC,sidenum);
    
    % Develop C-matrix from K-matrix (functions below)
    C = [];
    [C,~,~] = generateC(sel,rvar,NC,CA,Avar,E,C,sidenum);

    % Print C-matrix as output
    %disp('The C-matrix is: '); disp(C);
    
    % Calculate and print volume fraction (function below)
    volFrac = calcVF(NC,CA,rvar,Avar,sel,sidenum);
    %disp('The volume fraction is: '); disp(volFrac);
    
    % Plot design
    %plotDesign(NC,CA);
    
    % Plot nodal displacement (function below)
    %plotNDisp(NC,uBasket,sel);
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

% FUNCTION TO CALCULATE C-MATRIX
function [C,uBasket,FBasket] = generateC(sel,rvar,NC,CA,Avar,E,C,sidenum)
    % Initialize outputs
    uBasket = []; FBasket = [];
    
    % Iterate through once for each strain component
    for y = 1:1:3
    %   Define vectors to hold indexes for output forces
        Fi_x = []; Fi_y = []; Fi_xy = [];
    
    %   Define strain vector: [e11, e22, e12]'
        strainvec = [0;0;0];

    %   set that component equal to a dummy value (0.01 strain), 
    %       set all other values to zero
        strainvec(y) = 0.01; 
        strainvec(3) = strainvec(3)*2;

    %   use strain relations, BCs, and partitioned K-matrix to 
    %       solve for all unknowns
        e11 = strainvec(1); e22 = strainvec(2); e12 = strainvec(3);
        if (e11 ~= 0) || (e22 ~= 0) % when testing non-zero e11 or e22
            K = formK(NC,CA,Avar,E); % function for this below
            u_r = []; F_q = []; qvec = []; rvec = [];
            % Assigning Force/Displacement BCs for different nodes/DOFs
            for x = 1:1:size(NC,1) % looping through nodes by coordinate
                ND = NC./sel;
                % Separating conditions for exterior nodes
                if (ismember(ND(x,1),[0,1]) == true) || ...
                (ismember(ND(x,2),[0,1]) == true)
                    % Finding x-DOF
                    if ND(x,1) == 0
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((2*x)-1)];
                    elseif ND(x,1) == 1
                        % displacement in x = e11*sel
                        u_r = [u_r;(e11*sel)];
                        rvec = [rvec,((2*x)-1)];
                        Fi_x = [Fi_x,((2*x)-1)];
                    elseif (ND(x,2) == 0) && (e22 ~= 0)
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((2*x)-1)];
                    else
                        % x-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((2*x)-1)];
                    end
                    
                    % Finding y-DOF
                    if ND(x,2) == 0
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(2*x)];
                    elseif ND(x,2) == 1
                        % displacement in y = e22*sel
                        u_r = [u_r;(e22*sel)];
                        rvec = [rvec,(2*x)];
                        Fi_y = [Fi_y,(2*x)];
                        Fi_xy = [Fi_xy,((2*x)-1)];
                    elseif (ND(x,1) == 0) && (e11 ~= 0)
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((2*x)-1)];
                    else
                        % y-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,(2*x)];
                    end
                else % Condition for all interior nodes
                    % both x- and y-DOFs are force BCs
                    F_q = [F_q;0;0];
                    qvec = [qvec,((2*x)-1),(2*x)];
                end
            end
            qrvec = [qvec,rvec];
            newK = [K(qvec,:);K(rvec,:)];
            newK = [newK(:,qvec),newK(:,rvec)];
            K_qq = newK(1:length(qvec),1:length(qvec));
            K_rq = newK((length(qvec)+1):(2*size(NC,1)),1:length(qvec));
            K_qr = newK(1:length(qvec),(length(qvec)+1):(2*size(NC,1)));
            K_rr = newK((length(qvec)+1):(2*size(NC,1)),...
                   (length(qvec)+1):(2*size(NC,1)));
            u_q = K_qq\(F_q-(K_qr*u_r)); 
            F_r = (K_rq*u_q)+(K_rr*u_r);
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end
        else % when testing non-zero e12
            K = formK(NC,CA,Avar,E); % function for this below
            u_r = []; F_q = []; qvec = []; rvec = [];
            % Assigning Force/Displacement BCs for different nodes/DOFs
            for x = 1:1:size(NC,1) % looping through nodes by coordinate
                % Separating conditions for exterior nodes 
                ND = NC./sel;
                if (ismember(ND(x,1),[0,1]) == true) || ...
                (ismember(ND(x,2),[0,1]) == true)
                    % Finding x-DOF
                    if ND(x,1) == 0
                        % displacement in x is proportional to y-coordinate
                        % (due to nature of shear)
                        u_r = [u_r;(e12*sel*ND(x,2))];
                        rvec = [rvec,((2*x)-1)];
                    elseif ND(x,1) == 1
                        % displacement in x is proportional to y-coordinate
                        % (due to nature of shear)
                        u_r = [u_r;(e12*sel*ND(x,2))];
                        rvec = [rvec,((2*x)-1)];
                        Fi_x = [Fi_x,((2*x)-1)];
                    elseif ND(x,2) == 1
                        % displacement in x = e12*sel
                        u_r = [u_r;(e12*sel)];
                        rvec = [rvec,((2*x)-1)];
                    elseif ND(x,2) == 0
                        % displacement in x = 0
                        u_r = [u_r;0];
                        rvec = [rvec,((2*x)-1)];
                    else
                        % x-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,((2*x)-1)];
                    end
                    
                    % Finding y-DOF
                    if ND(x,2) == 0
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(2*x)];
                    elseif ND(x,2) == 1
                        % displacement in y = 0
                        u_r = [u_r;0];
                        rvec = [rvec,(2*x)];
                        Fi_y = [Fi_y,(2*x)];
                        Fi_xy = [Fi_xy,((2*x)-1)];
                    else
                        % y-direction DOF is a force BC
                        F_q = [F_q;0];
                        qvec = [qvec,(2*x)];
                    end
                else % Blanket condition for all interior nodes
                    % both x- and y-DOFs are force BCs
                    F_q = [F_q;0;0];
                    qvec = [qvec,((2*x)-1),(2*x)];
                end
            end
            qrvec = [qvec,rvec];
            newK = [K(qvec,:);K(rvec,:)];
            newK = [newK(:,qvec),newK(:,rvec)];
            K_qq = newK(1:length(qvec),1:length(qvec));
            K_rq = newK((length(qvec)+1):(2*size(NC,1)),1:length(qvec));
            K_qr = newK(1:length(qvec),(length(qvec)+1):(2*size(NC,1)));
            K_rr = newK((length(qvec)+1):(2*size(NC,1)),...
                   (length(qvec)+1):(2*size(NC,1)));
            u_q = K_qq\(F_q-(K_qr*u_r));
            F_r = (K_rq*u_q)+(K_rr*u_r);
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end 
        end
    
    %   Finding average side "thicknesses" due to differing element radii
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
        horizmean = mean(horizrads);
        vertmean = mean(vertrads);
    
    %   use system-wide F matrix and force relations to solve for 
    %       stress vector [s11,s22,s12]'
        F_x = 0; F_y = 0; F_xy = 0;
        for n = 1:1:size(Fi_xy,2)
            F_x = F_x + F(Fi_x(n));
            F_y = F_y + F(Fi_y(n));
            F_xy = F_xy + F(Fi_xy(n));
        end
        stressvec = [F_x/(sel*2*vertmean);F_y/(sel*2*horizmean);...
                     F_xy/(sel*2*horizmean)];
        if y == 1
            stresses = stressvec;
        end

    %   use strain and stress vectors to solve for the corresponding
    %       row of the C matrix
        Cdummy = stressvec/strainvec;
        C(:,y) = Cdummy(:,y);
        FBasket(:,y) = F;
        uBasket(:,y) = u;
    end
end

% FUNCTION TO FORM GLOBAL STRUCTURAL STIFFNESS MATRIX
function K = formK(NC,CA,Avar,E)
    % Forming Elemental Stiffness Matrices
    Kbasket = [];
    for i = 1:1:size(CA,1)
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        c=(x2-x1)/L; c2 = c^2;
        s=(y2-y1)/L; s2 = s^2;
        ktemp = [c2,   c*s,  -c2,  -c*s;  
                 c*s,   s2,  -c*s,  -s2;    
                 -c2, -c*s,  c2,    c*s; 
                 -c*s, -s2,  c*s,   s2];
        ke = ((Avar(i).*E)./L).*ktemp;
        Kbasket(:,:,i) = ke;
    end

    % Global-to-local-coordinate-system Coordination
    GlobToLoc=zeros(size(CA,1),4);
    for n=1:2  
        GN=CA(:,n); 
        for d=1:2
            GlobToLoc(:,(n-1)*2+d)=(GN-1)*2+d;
        end
    end

    % Forming Global Truss Stiffness Matrix
    K = zeros(2*size(NC,1));
    for e=1:size(CA,1) 
        ke = Kbasket(:,:,e);
        for lr = 1:4
            gr = GlobToLoc(e,lr); 
            for lc = 1:4
                gc = GlobToLoc(e,lc); 
                K(gr,gc) = K(gr,gc) + ke(lr,lc);
            end
        end
    end
end

% FUNCTION TO CALCULATE VOLUME FRACTION 
function volFrac = calcVF(NC,CA,rvar,Avar,sel,sidenum)
    totalTrussVol = 0;
    totl = sel;
    for i = 1:size(CA,1)
        % Finding element length from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        % Adding current element to total volume of trusses
        totalTrussVol = totalTrussVol + (L*Avar(i));
    end
    
    % Finding average side "thickness" due to differing element radii
    horizrads = [];
    for i = 1:1:size(CA,1)
        for j = 1:1:(sidenum-1)
            if ((CA(i,1) + (j*sidenum)) == CA(i,2)) && ...
                    ((NC(CA(i,1),2) == totl) || (NC(CA(i,1),2) == 0))
                singl = totl/(sidenum-1);
                x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
                y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
                L = sqrt(((x2-x1)^2)+((y2-y1)^2));
                for k = 1:1:(L/singl)
                    horizrads = [horizrads,rvar(i)];
                end
            elseif ((CA(i,1) - (j*sidenum)) == CA(i,2)) && ...
                    ((NC(CA(i,1),2) == totl) || (NC(CA(i,1),2) == 0))
                singl = totl/(sidenum-1);
                x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
                y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
                L = sqrt(((x2-x1)^2)+((y2-y1)^2));
                for k = 1:1:(L/singl)
                    horizrads = [horizrads,rvar(i)];
                end
            end
        end
    end
    vertrads = [];
    for i = 1:1:size(CA,1)
        for j = 1:1:(sidenum-1)
            if ((CA(i,1) + j) == CA(i,2)) && ...
                    ((NC(CA(i,1),1) == totl) || (NC(CA(i,1),1) == 0))
                singl = totl/(sidenum-1);
                x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
                y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
                L = sqrt(((x2-x1)^2)+((y2-y1)^2));
                for k = 1:1:(L/singl)
                    vertrads = [vertrads,rvar(i)];
                end
            elseif ((CA(i,1) - j) == CA(i,2)) && ...
                    ((NC(CA(i,1),1) == totl) || (NC(CA(i,1),1) == 0))
                singl = totl/(sidenum-1);
                x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
                y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
                L = sqrt(((x2-x1)^2)+((y2-y1)^2));
                for k = 1:1:(L/singl)
                    vertrads = [vertrads,rvar(i)];
                end    
            end
        end
    end
    thick = mean([mean(horizrads),mean(vertrads)]);
    
    % Calculating volume fraction (using a solid square with 2*(avg 
    %   thickness) as a baseline)
    volFrac = totalTrussVol/(2*thick*(totl^2)); 
end

% FUNCTION TO PLOT DESIGN
function plotDesign(NC,CA)
    for i = 1:1:size(CA,1)
        figure(1)
        plot([NC(CA(i,1),1),NC(CA(i,2),1)],[NC(CA(i,1),2),NC(CA(i,2),2)],'b-');
        hold on;
    end
end

% FUNCTION TO PLOT NODAL DISPLACEMENT
function plotNDisp(NB,uBasket,sel)
    NC = NB./sel;
    sf = 5; % scale factor (for magnifying displacement plotting)
    newNC = [];
    n = 1;
    for i = 1:1:size(uBasket,2)
        u = uBasket(:,i);
        nnc = [];
        n = 1;
        for j = 1:2:size(u)
            newrow = [NC(n,1)+(sf*u(j)),NC(n,2)+(sf*u(j+1))];
            nnc = [nnc;newrow];
            n = n+1;
        end
        newNC(:,:,i) = nnc;
    end
    subplot(2,2,1);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,1),newNC(:,2,1),'r*');
    %axis([-0.2 1.2 -0.2 1.2]);
    title('X-Direction Axial Strain (e11)');
    
    subplot(2,2,2);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,2),newNC(:,2,2),'g*');
    %axis([-0.2 1.2 -0.2 1.2]);
    title('Y-Direction Axial Strain (e22)');
    
    subplot(2,2,3);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,3),newNC(:,2,3),'k*');
    %axis([-0.2 1.2 -0.2 1.2]);
    title('Shear Strain (e12)');
end

