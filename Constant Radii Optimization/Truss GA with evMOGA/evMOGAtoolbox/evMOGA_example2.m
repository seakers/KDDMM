%% evMOGA example 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Minimal algorithm parameters set (problem characteristics)
clear eMOGA
eMOGA.objfun='mop3';           % m-function name for objectives computation
eMOGA.objfun_dim=2;            % Objective space dimension
eMOGA.searchspaceUB=[pi pi];   % Search space upper bound
eMOGA.searchspaceLB=[-pi -pi]; % Search space lower bound
eMOGA.Nind_P= 500;             % Individuals for the P population
eMOGA.Generations= 100;        % Number of generations
eMOGA.n_div= 200;              % Number of division for each dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm execution
[pfront,pset,eMOGA]=evMOGA(eMOGA);
