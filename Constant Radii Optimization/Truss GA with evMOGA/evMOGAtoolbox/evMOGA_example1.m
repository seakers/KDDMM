%% evMOGA example 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Minimal algorithm parameters set (problem characteristics)
clear eMOGA
eMOGA.objfun='mop3';            % m-function name for objectives computation
eMOGA.objfun_dim=2;             % Objective space dimension
eMOGA.searchspaceUB=[pi pi];    % Search space upper bound
eMOGA.searchspaceLB=[-pi -pi];  % Search space lower bound

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm execution
[pfront,pset,eMOGA]=evMOGA(eMOGA);
