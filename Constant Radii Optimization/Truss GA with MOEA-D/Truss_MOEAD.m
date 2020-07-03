%% Truss MOEA-D
clear;  
close all;
clc;

%% Problem Constants
E = 10000; % Young's Modulus for polymeric material (example: 10000 Pa)
sel = 0.05; % Unit square side length (NOT individual truss length) (example: 5 cm)
r = 50*(10^-6); % Radius for cross-sectional area of (assumed circular) truss members (example: 50 micrometers)
A = pi*(r^2); % Cross-sectional area of truss member
% Nodal Coordinate Vector (Standard 3x3, 2D Grid) below (each row represents a node, first column is x-coordinates, second column is y-coordinates):
NC = sel.*[0,0;0,0.5;0,1;0.5,0;0.5,0.5;0.5,1;1,0;1,0.5;1,1]; 
CA_all = [1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 2,3; 2,4; 2,5; 2,6; 2,7; 2,8; 2,9; 3,4; 3,5; 3,6; 3,7; 3,8; 3,9; 4,5; 4,6; 4,7; 4,8; 4,9; 5,6; 5,7; 5,8; 5,9; 6,7; 6,8; 6,9; 7,8; 7,9; 8,9];
c_ratio = 1; % Ratio of C22/C11 to target
global CA_infeas;
CA_infeas = {};
global inf_con;
inf_con = 0;

%% MOEA-D Problem Definition
CostFunction=@(x) multiobjective_MOEAD(x,CA_all,NC,A,E,sel,r,c_ratio,CA_infeas,inf_con);  % Cost Function
nVar=36;             % Number of Decision Variables
VarSize=[nVar 1];   % Decision Variables Matrix Size
VarMin = 0;         % Decision Variables Lower Bound
VarMax = 1;         % Decision Variables Upper Bound
nObj=2;

%% MOEA/D Settings

MaxIt=100;  % Maximum Number of Iterations

nPop=500;    % Population Size (Number of Sub-Problems)

nArchive=50;

T=max(ceil(0.15*nPop),2);    % Number of Neighbors
T=min(max(T,2),15);

crossover_params.gamma=0.5;
crossover_params.VarMin=VarMin;
crossover_params.VarMax=VarMax;

%% Initialization

% Create Sub-problems
sp=CreateSubProblems(nObj,nPop,T);

% Empty Individual
empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.g=[];
empty_individual.IsDominated=[];

% Initialize Goal Point
%z=inf(nObj,1);
z=zeros(nObj,1);

% Create Initial Population
pop=repmat(empty_individual,nPop,1);
for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Cost=CostFunction(pop(i).Position);
    z=min(z,pop(i).Cost);
end

for i=1:nPop
    pop(i).g=DecomposedCost(pop(i),z,sp(i).lambda);
end

% Determine Population Domination Status
pop=DetermineDomination(pop);

% Initialize Estimated Pareto Front
EP=pop(~[pop.IsDominated]);

%% Main Loop

for it=1:MaxIt
    for i=1:nPop
        
        % Reproduction (Crossover)
        K=randsample(T,2);
        
        j1=sp(i).Neighbors(K(1));
        p1=pop(j1);
        
        j2=sp(i).Neighbors(K(2));
        p2=pop(j2);
        
        y=empty_individual;
        y.Position=Crossover(p1.Position,p2.Position,crossover_params);
        
        y.Cost=CostFunction(y.Position);
        
        z=min(z,y.Cost);
        
        for j=sp(i).Neighbors
            y.g=DecomposedCost(y,z,sp(j).lambda);
            if y.g<=pop(j).g
                pop(j)=y;
            end
        end
        
    end
    
    % Determine Population Domination Status
	pop=DetermineDomination(pop);
    
    ndpop=pop(~[pop.IsDominated]);
    
    EP=[EP
        ndpop]; %#ok
    
    EP=DetermineDomination(EP);
    EP=EP(~[EP.IsDominated]);
    
    if numel(EP)>nArchive
        Extra=numel(EP)-nArchive;
        ToBeDeleted=randsample(numel(EP),Extra);
        EP(ToBeDeleted)=[];
    end
    
    % Plot EP
    figure(1);
    PlotCosts(EP);
    pause(0.01);

    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Number of Pareto Solutions = ' num2str(numel(EP))]);
    
end

%% Reults

disp(' ');

EPC=[EP.Cost];
for j=1:nObj
    
    disp(['Objective #' num2str(j) ':']);
    disp(['      Min = ' num2str(min(EPC(j,:)))]);
    disp(['      Max = ' num2str(max(EPC(j,:)))]);
    disp(['    Range = ' num2str(max(EPC(j,:))-min(EPC(j,:)))]);
    disp(['    St.D. = ' num2str(std(EPC(j,:)))]);
    disp(['     Mean = ' num2str(mean(EPC(j,:)))]);
    disp(' ');
    
end

%% Compute and plot feasibility and stability scores for final population
archive = [EP.Position];
feas_scores = zeros(size(archive,2),1);
stab_scores = zeros(size(archive,2),1);
sidenum = 3;
for i = 1:size(archive,2)
    x_curr = archive(:,i).';
    
    % Converting design vector to CA matrix
    x_bin = x_curr>0.5;
    CA_des = CA_all(x_bin~=0,:);
    
    %CA_des = CA_all(x_curr~=0,:);
    
    %x_des = [0, 1, 0, 0, x_curr(1:2), 1, x_curr(3:4), 0, x_curr(5:18), ...
    %x_curr(19:end), 0, 1, 0]; % for forcing truss elements in 1,3; 1,7 and
    % 7;9
    %CA_des = CA_all(x_des~=0,:);
    
    %x_des = [x_curr(1), 0, x_curr(2:3), 0, 0, 0, 0, x_curr(4:7), 0, 0, 0, 0, ... 
        %x_curr(8:9), 0, 0, 0, x_curr(10), 0, x_curr(11:12), 0, x_curr(13:16), ...
        %0, x_curr(17:19), 0, x_curr(20)]; % only adjacent node connections allowed
    %CA_des = CA_all(x_des~=0,:);
    
    % Computing feasibility and stability scores
    feas_scores(i) = feasibility_checker_nonbinary(NC,CA_des);
    stab_scores(i) = stabilityTester_2D_updated(sidenum,CA_des,NC);
end

% Plotting
figure()
plot(feas_scores, stab_scores, 'b*')
xlabel('Feasibility scores')
ylabel('Stability scores')
title('Constraints comparison - Method 5')

%% Plotting Pareto Front (SEAK code)
%%% Finding only feasible designs
objectives = [EP.Cost];
feasibility = false(size(archive,2),1);
feasibility_nonbinary = zeros(size(archive,2),1);
for i = 1:size(archive,2)
    x_curr = archive(:,i).';
    [feasibility_nonbinary(i), feasibility(i)] = feasibility_checker_boolean(x_curr, NC, CA_all); 
end
x_feas = archive(:,feasibility).';
f_feas = objectives(:,feasibility).';
%plot_pareto_seak(f_feas)

% Plotting the true feasible objectives
f_true_feas = zeros(size(x_feas,2),2);
lambda = 100;
for i = 1:size(x_feas,1)
    f_true_feas(i,:) = f_feas(i,:) + lambda*2*ones(1,2);  
end
plot_pareto_seak(f_true_feas,2)

%%% Plotting all designs
%plot_pareto_seak(pfront)

%% Convert design vector to binary array (not needed for bitstring design vector)
x_bin_feas = false(size(x_feas,1), size(CA_all,1));
for i = 1:size(x_feas,1)
    x_feas_i = x_feas(i,:);
    x_bin_feas(i,:) = x_feas_i>0.5; 
end

%% Visualize trusses
x_feas_unique = unique(x_feas,'rows');
for i = 1:size(x_feas_unique,1)
    visualize_truss_fromx_3x3(NC, CA_all, x_feas_unique(i,:))
end