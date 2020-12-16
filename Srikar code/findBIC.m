% Generating Data from Roshan's GA runs to single text file to find BIC for
% best fit determination

% Initialization
clc;   
close all; 
clear;

UC = [1,2,3,4,6,8,10]; UC2 = ([1,2,3,4,6,8,10].^2);
C11basket = []; C12basket = []; C44basket = [];

% Input CAs from GA runs
CABasket = extractDesigns();

% Restricting CA's considered by size
CA_alldesigns = CABasket;
row = cellfun(@(x) size(x,1),CA_alldesigns);
CABasket = CA_alldesigns(row>12);

% Loop through each design to solve for, plot C-matrix values
for cac = 1:1:size(CABasket,1)
    tinyCA = cell2mat(CABasket(cac));
    %disp(tinyCA);
    C11 = []; C12 = []; C44 = [];
    
    % Loop through each # of unit cells
    for nuc = 1:1:size(UC,2)
        % Input Constants
        nucFac = UC(nuc); % factor to indicate how many unit cells are in  
              % this model (nucFac^2 gives the actual number of unit cells)
        sel = 0.05; % unit cube side length of 5 cm, truss length of 2.5 cm
        r = 50*(10^-6); % Radius of 50 micrometers for cross-sectional 
                         %  area of (assumed circular) truss members
        E = 10000; % Young's Modulus of 10000 Pa (polymeric material)

        % Calculate C, volFrac for each # of unit cells case
        
        [C,volFrac] = calcCTruss(nucFac,sel,r,E,tinyCA);

        % Record outputs for plotting/reference
        C11 = [C11,C(1,1)];
        %C12 = [C12,C(1,2)];
        %C44 = [C44,C(3,3)];
    end
    C11basket = [C11basket;C11]; 
    %C12basket = [C12basket,C12]; 
    %C44basket = [C44basket,C44]; 
end   

% Calculate information criterion for each design
aicExp1Basket = []; bicExp1Basket = [];
aicExp2Basket = []; bicExp2Basket = [];

% Loop through results of each design
for k = 1:1:size(C11basket,1)
    Cvals = C11basket(k,:);
    
    % Find fits
    fexp1 = fit(UC',Cvals','exp1'); f1 = coeffvalues(fexp1);
    x1 = f1(1).*exp(f1(2).*UC');
    fexp2 = fit(UC',Cvals','exp2'); f2 = coeffvalues(fexp2);
    x2 = f2(1).*exp(f2(2).*UC')+f2(3).*exp(f2(4).*UC');
    
    % Find loglikelihoods
    mu = 1;
    %sigma = 1;
    pd = makedist('Exponential','mu',mu);
    logL = [sum(log(pdf(pd,x1))),sum(log(pdf(pd,x2)))];

    % Input number of models, number of observations
    numParam = length(UC).*ones(1,2);
    numObs = size(CABasket,1);

    % Find Information Criteria
    [aic,bic] = aicbic(logL,numParam,numObs);
    aicExp1Basket(k) = aic(1); bicExp1Basket(k) = bic(1);
    aicExp2Basket(k) = aic(2); bicExp2Basket(k) = bic(2);
end

% Find average information criterion
avgaicExp1 = mean(aicExp1Basket); avgbicExp1 = mean(bicExp1Basket);
avgaicExp2 = mean(aicExp2Basket); avgbicExp2 = mean(bicExp2Basket);

%----------%
% FUNCTION TO EXTRACT ALL DESIGNS FROM SPREADSHEETS
function CA_alldesigns = extractDesigns()
    filename1 = "EpsilonMOEA_emoea_";
    filename2 = "fscon0__fibre.csv";

    n_runs = 30;

    designs_feas_and_stab_allruns = [];

    for i = 1:n_runs
        full_filepath = strcat(filename1,num2str(i-1),filename2);
        data_table = readtable(full_filepath,'Format','%s%f%f%f%f',...
            'HeaderLines',1);

        %%%% store retrieved data into different variables
        %%%% csv_data includes: [Pen. Obj. 1, Pen.Obj. 2,Feasibility Score,
        %%%% Stablity Score]
        pop_size =  size(data_table,1);
        csv_data = zeros(pop_size,4);
        designs = strings(pop_size);
        csv_data = data_table(:,2:end);
        designs = data_table(:,1);

        csv_data_array = table2array(csv_data);
        designs_array = table2array(designs);

        feas_array = csv_data_array(:,3);
        stab_array = csv_data_array(:,4);

        designs_array_feas_and_stab = designs_array((feas_array == 1) & ...
            (stab_array == 1));
        designs_feas_and_stab_allruns = [designs_feas_and_stab_allruns;...
            designs_array_feas_and_stab];
    end

    % Remove non-unique designs 
    designs_allruns_unique = unique(designs_feas_and_stab_allruns);

    % Convert to Connectivity Array 
    % Boolean string representation of design is converted to CA matrix and
    % stored in a structure
    CA_all = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 2,3; 2,4; 2,5;2,6;...
        2,7; 2,8; 2,9; 3,4; 3,5; 3,6; 3,7; 3,8; 3,9; 4,5; 4,6; 4,7; 4,8;...
        4,9; 5,6; 5,7; 5,8; 5,9; 6,7; 6,8; 6,9; 7,8; 7,9; 8,9];
    CA_allruns_struct = struct; 

    for i = 1:size(designs_allruns_unique)
        current_des_boolstring = designs_allruns_unique{i};
        design_bool = zeros(size(current_des_boolstring,2),1);
        for j = 1:size(current_des_boolstring,2)
            design_bool(j) = str2num(current_des_boolstring(j));
        end

        CA_des = CA_all(design_bool~=0,:);
        current_field = strcat('design',num2str(i));
        CA_allruns_struct.(current_field) = CA_des;
    end
    
    CA_alldesigns = struct2cell(CA_allruns_struct);
end

% FUNCTION TO CALCULATE FULL CA, THEN C, FOR GIVEN # OF UNIT CELLS
function [C,volFrac] = calcCTruss(nucFac,sel,r,E,tinyCA)
    % Calculated Inputs
    sidenum = (2*nucFac) + 1; % Number of nodes on each size of square grid
    A = pi*(r^2); % Constant cross-sectional area of truss members

    % Generate vector with nodal coordinates
    NC = generateNC(sel,sidenum);

    % Generate Connectivity Array for given # of UC's
    CA = generateCA(sidenum,tinyCA);

    % Develop C-matrix from K-matrix (functions below)
    C = [];
    [C,~,~] = generateC(sel,r,NC,CA,A,E,C);
    C = C./nucFac;

    % Calculate volume fraction (function below)
    volFrac = calcVF(NC,CA,r,sel);
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

% FUNCTION TO GENERATE CONNECTIVITY ARRAY FOR 3x3 UC'S, nucFac x nucFac 
function CA = generateCA(sidenum,tinyCA)
    CA = [];

    % Outside frame
    for i = 1:1:(sidenum-1)
        CA = [CA;i,(i+1)];
    end
    for i = sidenum:sidenum:((sidenum^2)-sidenum)
        CA = [CA;i,(i+sidenum)];
    end
    for i = ((sidenum^2)-sidenum+1):1:((sidenum^2)-1)
        CA = [CA;i,(i+1)];
    end
    for i = 1:sidenum:((sidenum^2)-(2*sidenum)+1)
        CA = [CA;i,(i+sidenum)];
    end
    
    % Identify central node for all unit cells
    cnvec = [];
    for j = (sidenum+2):(2*sidenum):((sidenum^2)-(2*sidenum)+2)
        cnvec = [cnvec;j];
        for i = 2:2:(sidenum-2)
            cnvec = [cnvec;(j+i)];
        end
    end
    
    % Populate CA using template of tinyCA and locations of cn's
    for w = 1:1:size(cnvec,1)
        cn = cnvec(w);
        for q = 1:1:size(tinyCA,1)
            row = [];
            for r = 1:1:2
                if tinyCA(q,r) == 1
                    row = [row,(cn-sidenum-1)];
                elseif tinyCA(q,r) == 2
                    row = [row,(cn-sidenum)];
                elseif tinyCA(q,r) == 3
                    row = [row,(cn-sidenum+1)];
                elseif tinyCA(q,r) == 4
                    row = [row,(cn-1)];
                elseif tinyCA(q,r) == 5
                    row = [row,cn];
                elseif tinyCA(q,r) == 6
                    row = [row,(cn+1)];
                elseif tinyCA(q,r) == 7
                    row = [row,(cn+sidenum-1)];
                elseif tinyCA(q,r) == 8
                    row = [row,(cn+sidenum)];
                elseif tinyCA(q,r) == 9
                    row = [row,(cn+sidenum+1)];
                end
            end
            CA = [CA;row];
        end
    end
    
    % Eliminate duplicate members
    CA = unique(CA,'rows');
end

% FUNCTION TO CALCULATE C-MATRIX
function [C,uBasket,FBasket] = generateC(sel,r,NC,CA,A,E,C)
    % Initialize outputs
    uBasket = []; FBasket = [];

    % Iterate through once for each strain component
    for y = 1:1:3
      % Define vectors to hold indexes for output forces
        Fi_x = []; Fi_y = []; Fi_xy = [];

      % Define strain vector: [e11, e22, e12]'
        strainvec = [0;0;0];

      % set that component equal to a dummy value (0.01 strain), 
      %     set all other values to zero
        strainvec(y) = 0.01; 
        strainvec(3) = strainvec(3)*2;

      % use strain relations, BCs, and partitioned K-matrix to 
      %     solve for all unknowns
        e11 = strainvec(1); e22 = strainvec(2); e12 = strainvec(3);
        if (e11 ~= 0) || (e22 ~= 0) % when testing non-zero e11 or e22
            K = formK(NC,CA,A,E); % function for this below
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
            K = formK(NC,CA,A,E); % function for this below
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

      % use system-wide F matrix and force relations to solve for 
      %     stress vector [s11,s22,s12]'
        F_x = 0; F_y = 0; F_xy = 0;
        for n = 1:1:size(Fi_xy,2)
            F_x = F_x + F(Fi_x(n));
            F_y = F_y + F(Fi_y(n));
            F_xy = F_xy + F(Fi_xy(n));
        end
        stressvec = (1/(sel*2*r)).*[F_x;F_y;F_xy];

      % use strain and stress vectors to solve for the corresponding
      %     row of the C matrix
        Cdummy = stressvec/strainvec;
        C(:,y) = Cdummy(:,y);
        FBasket(:,y) = F;
        uBasket(:,y) = u;
    end
end

% FUNCTION TO FORM GLOBAL STRUCTURAL STIFFNESS MATRIX
function K = formK(NC,CA,A,E)
    % Forming Elemental Stiffness Matrices
    Kbasket = [];
    for i = 1:size(CA,1)
        x1 = NC(CA(i,1),1); x2 = NC(CA(i,2),1);
        y1 = NC(CA(i,1),2); y2 = NC(CA(i,2),2);
        L = sqrt(((x2-x1)^2)+((y2-y1)^2));
        c=(x2-x1)/L; c2 = c^2;
        s=(y2-y1)/L; s2 = s^2;
        ktemp = [c2,   c*s,  -c2,  -c*s;  
                 c*s,   s2,  -c*s,  -s2;    
                 -c2, -c*s,  c2,    c*s; 
                 -c*s, -s2,  c*s,   s2];
        ke = ((A.*E)./L).*ktemp;
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
