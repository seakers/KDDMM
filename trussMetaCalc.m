% Function for calculating metamaterial properties for a 2D truss
% All lengths are in [m], all stresses and moduli are in [Pa] 
% The nodes on the 3x3 grid are numbered as such:
%       3   6   9
%       2   5   8  
%       1   4   7
% Each node has two degrees of freedom (DOF), x and y.  The (node number*2) 
%    gives the number for that node's yDOF, and (yDOF - 1) is the xDOF
% We assume that all edge elements are present

function [C,stab_pre,symm_post,volFrac] = trussMetaCalc(sel,r,E,CA)

%-% Default Inputs (leave commented out)
%sel = 0.05; % Unit square side length (NOT individual truss length)
             %  (example: 5 cm)
%r = 50*(10^-6); % Radius for cross-sectional area of (assumed circular)
                 %   truss members (example: 50 micrometers)
%E = 10000; % Young's Modulus for polymeric material (example: 10000 Pa)
% Connectivity array (used to indicate position of elements): each row
%   contains the two nodes that an element straddles
%CA = [1,2;2,3;1,4;1,5;2,5;3,5;3,6;4,5;5,6;4,7;5,7;5,8;5,9;6,9;7,8;8,9];


%-% Fixed Inputs
A = pi*(r^2); % Cross-sectional area of truss member
% Nodal Coordinate Vector (Standard 3x3, 2D Grid) below (each row
%   represents a node, first column is x-coordinates, second column is
%   y-coordinates):
NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1]; 

%-% Check stability pre-FEA (function below)
% N is the nodal connectivity vector: each index of N represents a node, 
%   and its value represents the # of elements connected to it
[N,stab_pre] = stabilityTester(CA);
if stab_pre == 0
    disp('This truss design is not mechanically stable');
    return
end

%-% Develop C-matrix from K-matrix (functions below)
C = []; % units are [Pa]
[C,uBasket,FBasket] = generateC(sel,r,NC,CA,A,E,C);

%-% Calculate volume fraction (function below)
volFrac = calcVF(NC,CA,r,sel);

%-% Plot nodal displacement (function below)
%plotNDisp(NC,uBasket);

%-% Print C-matrix as output
%disp('The C-matrix is: '); disp(C);

%-% Post-FEA C-symmetry check (function below)
symm_post = symmCCheck(C);
if symm_post == 0
    disp('This truss design is not eligible due to a lack of homogenizability');
end

end


%----------%
% FUNCTION TO TEST TRUSS STABILITY
function [N,stability] = stabilityTester(CA)
    stability = [];

    % Add up nodal connectivity vector N: each index of N represents a
    %   node, and its value represents the # of elements connected to it
    [N,~] = histcounts(CA,9);
    
    % Determine stability based on values in N
    if N(5)>=3 % Case 1: stability is tied to central node connections
        if (N(1)>=3)||(N(7)>=3)||(N(3)>=3)||(N(9)>=3)
            if (N(1)>=3)&&(N(3)>=3)
                stability = true;
            elseif (N(7)>=3)&&(N(3)>=3)
                stability = true;
            elseif (N(9)>=3)&&(N(3)>=3)
                stability = true;
            elseif (N(1)>=3)&&(N(7)>=3)
                stability = true;
            elseif (N(1)>=3)&&(N(9)>=3)
                stability = true;
            elseif (N(7)>=3)&&(N(9)>=3)
                stability = true;
            elseif (N(2)>=3)||(N(4)>=3)||(N(6)>=3)||(N(8)>=3)
                if (N(2)>=3)&&(N(4)>=3)
                    stability = true;
                elseif (N(6)>=3)&&(N(4)>=3)
                    stability = true;
                elseif (N(8)>=3)&&(N(4)>=3)
                    stability = true;
                elseif (N(8)>=3)&&(N(2)>=3)
                    stability = true;
                elseif (N(8)>=3)&&(N(6)>=3)
                    stability = true;
                elseif (N(6)>=3)&&(N(2)>=3)
                    stability = true;
                else
                    stability = false;
                end
            else
                stability = false;
            end  
        end
    % Case 2: stability isn't tied to central node connections
    elseif (N(2)>=4)||(N(4)>=4)||(N(6)>=4)||(N(8)>=4)
        if (N(2)>=4)&&(N(4)>=4)&&(N(6)>=4)
            stability = true;
        elseif (N(4)>=4)&&(N(6)>=4)&&(N(8)>=4)
            stability = true;
        elseif (N(2)>=4)&&(N(6)>=4)&&(N(8)>=4)
            stability = true;
        elseif (N(2)>=4)&&(N(4)>=4)&&(N(8)>=4)
            stability = true;
        else
            stability = false;
        end
    else
        stability = false;
    end
end

% FUNCTION TO FORM GLOBAL STRUCTURAL STIFFNESS MATRIX
function K = formK(NC,CA,A,E)
    % Forming Elemental Stiffness Matrices
    Kbasket = [];
    for i = 1:size(CA,1)
        % Finding element length from nodal coordinates
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        % Calculating directional cosines
        c=(x2-x1)/L; c2 = c^2;
        s=(y2-y1)/L; s2 = s^2;
        % Calculating elemental stiffness matrix
        ktemp = [c2,   c*s,  -c2,  -c*s;  
                 c*s,   s2,  -c*s,  -s2;    
                 -c2, -c*s,  c2,    c*s; 
                 -c*s, -s2,  c*s,   s2];
        ke = ((A.*E)./L).*ktemp;
        % Storing the elemental K matrices separately
        Kbasket(:,:,i) = ke;
    end

    % Global-to-local-coordinate-system coordination:
    %    used to identify the position of a given local node 
    %    for a given elemental K matrix in the global truss, 
    %    based on the relevant node number
    GlobToLoc=zeros(size(CA,1),4);
    for n=1:2  
        GN=CA(:,n); 
        for d=1:2
            GlobToLoc(:,(n-1)*2+d)=(GN-1)*2+d;
        end
    end

    % Forming Global Truss Stiffness Matrix:
    %    done by adding the local K matrices at each 
    %    degree of freedom, based on the global-to-
    %    local coordination
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

% FUNCTION TO CALCULATE C-MATRIX
function [C,uBasket,FBasket] = generateC(sel,r,NC,CA,A,E,C)
    % Initialize outputs
    uBasket = []; FBasket = [];
    
    % Iterate through once for each strain component
    for y = 1:1:3
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
            % Vector used to rearrange degrees of freedom
            qrvec = [4,7,9,10,11,16,1,2,3,5,6,8,12,13,14,15,17,18];
            % Rearranging K for rearranged degrees of freedom
            K = formK(NC,CA,A,E);
            newK = [K([4,7,9,10,11,16],:);K([1:3,5,6,8,12:15,17,18],:)];
            newK = [newK(:,[4,7,9,10,11,16]),...
                    newK(:,[1:3,5,6,8,12:15,17,18])];
            K = newK;
            % Defining known displacement, force boundary conditions
            u_r = [0;0;0;0;e22*sel;0;e22*sel;e11*sel;...
                   0;e11*sel;e11*sel;e22*sel];
            F_q = [0;0;0;0;0;0];
            % Partitioning K
            K_qq = K(1:6,1:6);
            K_rq = K(7:18,1:6);
            K_qr = K(1:6,7:18);
            K_rr = K(7:18,7:18);
            % Solving for unknown forces, displacements
            u_q = K_qq\(F_q-(K_qr*u_r)); 
            F_r = (K_rq*u_q)+(K_rr*u_r);
            % Assembling F,u vector in original order
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end
        else % when testing non-zero e12
            % Vector used to rearrange degrees of freedom
            qrvec = [4,9,10,16,1,2,3,5,6,7,8,11,12,13,14,15,17,18];
            % Rearranging K for rearranged degrees of freedom
            K = formK(NC,CA,A,E);            
            newK = [K([4,9,10,16],:);K([1:3,5:8,11:15,17,18],:)];
            newK = [newK(:,[4,9,10,16]),newK(:,[1:3,5:8,11:15,17,18])];
            K = newK;
            % Defining known displacement, force boundary conditions
            u_r = [0;0;0.5*e12*sel;e12*sel;0;0;0;e12*sel;0;0;...
                   0;0.5*e12*sel;e12*sel;0];
            F_q = [0;0;0;0];
            % Partitioning K
            K_qq = K(1:4,1:4);
            K_rq = K(5:18,1:4);
            K_qr = K(1:4,5:18);
            K_rr = K(5:18,5:18);
            % Solving for unknown forces, displacements
            u_q = K_qq\(F_q-(K_qr*u_r));
            F_r = (K_rq*u_q)+(K_rr*u_r);
            % Assembling F,u vector in original order
            altu = [u_q;u_r]; altF = [F_q;F_r];
            F = zeros(length(altF),1); u = zeros(length(altu),1);
            for x = 1:1:length(qrvec)
                F(qrvec(x)) = altF(x);
                u(qrvec(x)) = altu(x);
            end 
        end
        
    %   use F vector and geometry to solve for stress vector [s11,s22,s12]'
        F_x = F(13)+F(15)+F(17);
        F_y = F(6)+F(12)+F(18);
        F_xy = F(5)+F(11)+F(17);
        stressvec = (1/(sel*2*r)).*[F_x;F_y;F_xy];

    %   use strain and stress vectors to solve for the corresponding
    %       row of the C matrix
        Cdummy = stressvec/strainvec;
        C(:,y) = Cdummy(:,y);
        FBasket(:,y) = F;
        uBasket(:,y) = u;
    end
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


% FUNCTION TO PLOT NODAL DISPLACEMENT
function plotNDisp(NC,uBasket)
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
    axis([-0.01 0.06 -0.01 0.06]);
    title('X-Direction Axial Strain (e11)');
    
    subplot(2,2,2);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,2),newNC(:,2,2),'g*');
    axis([-0.01 0.06 -0.01 0.06]);
    title('Y-Direction Axial Strain (e22)');
    
    subplot(2,2,3);
    plot(NC(:,1),NC(:,2),'b.',newNC(:,1,3),newNC(:,2,3),'k*');
    axis([-0.01 0.06 -0.01 0.06]);
    title('Shear Strain (e12)');
end

% FUNCTION TO TEST C-MATRIX SYMMETRY
function symm_post = symmCCheck(C)
    if C(1,2) == C(2,1)
        if C(3,1) == C(1,3)
            if C(2,3) == C(3,2)
                symm_post = true;
            else
                symm_post = false;
            end
        else
            symm_post = false;
        end
    else
        symm_post = false;
    end
end
