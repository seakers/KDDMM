function eMOGA=evMOGAresults(eMOGA)
%   eMOGA=evMOGAresult(eMOGA)
%     User tasks to be evaluated at the end of the algorithm.
%     For instance to plot final results, final analysis, etc.

disp('------------------------------------------------')
disp('######   RESULT   #########')
disp(['Nind_A=' num2str(eMOGA.Nind_A)])
disp(['Objectives min values: ' num2str(eMOGA.min_f)])
disp('------------------------------------------------')

% Plotting results
load mop3aux          % load true pareto front
ParetoFront=eMOGA.coste_A(1:eMOGA.Nind_A,:);
% Plot true pareto front
plot(pfront(:,1),pfront(:,2),'.r') 
% Pareto front obtained by evMOGA
hold on,plot(ParetoFront(:,1),ParetoFront(:,2),'o')