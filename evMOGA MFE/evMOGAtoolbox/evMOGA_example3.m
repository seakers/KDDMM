%% evMOGA example 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Minimal algorithm parameters set (problem characteristics)
clear eMOGA
eMOGA.objfun='mop3_p';         % m-function name for objectives computation
eMOGA.searchspaceUB=[pi pi];   % Search space upper bound
eMOGA.searchspaceLB=[-pi -pi]; % Search space lower bound
eMOGA.objfun_dim=2;            % Objective space dimension
% additional parameter of the objective function
% in this case just to waist time 
p.n=1000000;
p.p=50000;
eMOGA.param=p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm execution
[pfront,pset,eMOGA]=evMOGA(eMOGA);
