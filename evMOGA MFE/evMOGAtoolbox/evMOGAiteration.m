function eMOGA=evMOGAiteration(eMOGA)
%   evMOGA=evMOGAiteration(evMOGA)
%     User tasks to be evaluated at the end of each iteration.
%     For instance to plot partial results, made a backup, etc.

disp(['######### generation ' num2str(eMOGA.gen_counter)])
disp(['Nind_A=' num2str(eMOGA.Nind_A)])
disp(['Objectives min values: ' num2str(eMOGA.min_f)])

