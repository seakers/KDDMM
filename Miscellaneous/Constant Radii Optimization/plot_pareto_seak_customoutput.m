function [pareto_bool_gen] = plot_pareto_seak_customoutput(f_vals, n_viable)
% This function uses paretofront.m, paretofront.c and paretofront.mexw64 to 
% determine the designs in the pareto front and plots them. This function 
% is customized for the custom output case wherein the designs in 
% different generations are to be plotted (add figure
% statement in main script) (ONLY FOR 2 OBJECTIVE CASE)
% Inputs: Objective values for all evaluated designs (n_pop x 2 x n_gen)
%         Viable count for each generation n_viable (n_gen x 1)
% Output: Boolean non-dominance array pareto_bool_gen (n_pop x floor(n_gen/10)+2)

gen_div = 15; % number of generations between plotting 

n_gen = size(n_viable,1);
n_pop = size(f_vals,1);
pareto_bool_gen = false(n_pop,floor(n_gen/gen_div)+2);

color = linspace(1,floor(n_gen/gen_div)+2,floor(n_gen/gen_div)+2);
current_gen_array = f_vals(:,:,1);
viable_count = n_viable(1,1);
f_current_gen = current_gen_array(1:viable_count,:);

figure
pareto_bool_1 = paretofront(f_current_gen);
pareto_bool_gen(:,1) = [pareto_bool_1; false(n_pop-viable_count,1)];
f_pareto_1 = f_current_gen(pareto_bool_1==1,:);
f_pareto_true_1 = [f_pareto_1(:,1), -f_pareto_1(:,2)]; % objective 1 is to be
% minimized while objective 2 is to be maximized
color1_all = ones(size(f_current_gen,1),1).*color(1); 
color1 = ones(size(f_pareto_true_1,1),1).*color(1);
scatter(f_current_gen(:,1),-f_current_gen(:,2),color1_all)
hold on
scatter(f_pareto_true_1(:,1),f_pareto_true_1(:,2),20,color1,'filled')
legendarg = strings((floor(n_gen/gen_div)+2)*2,1);
legendarg(1,1) = 'Generation 0'; 
legendarg(2,1) = 'Generation 0 - non-dominant';

hold on
for i = 1:floor(n_gen/gen_div)
    current_gen_array = f_vals(:,:,i*gen_div);
    viable_count = n_viable(i*gen_div,1);
    f_current_gen = current_gen_array(1:viable_count,:);
    
    pareto_bool_i = paretofront(f_current_gen);
    pareto_bool_gen(:,i+1) = [pareto_bool_i; false(n_pop-viable_count,1)];
    f_pareto_i = f_current_gen(pareto_bool_i==1,:);
    f_pareto_true_i = [f_pareto_i(:,1), -f_pareto_i(:,2)]; % objective 1 is to be
    % minimized while objective 2 is to be maximized  
    color_i_all = ones(size(f_current_gen,1),1).*color(i+1); 
    color_i = ones(size(f_pareto_true_i,1),1).*color(i+1);
    scatter(f_current_gen(:,1),-f_current_gen(:,2),20,color_i_all)
    hold on
    scatter(f_pareto_true_i(:,1),f_pareto_true_i(:,2),20,color_i,'filled')
    legendarg(2*i+1,1) = strcat('Generation ',string(i*gen_div)); 
    legendarg(2*i+2,1) = strcat('Generation ',string(i*gen_div),' - non-dominant');

    %current_pareto_bool = plot_pareto_seak(f_current_gen,2);
    hold on
end

current_gen_array = f_vals(:,:,end);
viable_count = n_viable(end,1);
f_current_gen = current_gen_array(1:viable_count,:);

pareto_bool_end = paretofront(f_current_gen);
pareto_bool_gen(:,end) = [pareto_bool_end; false(n_pop-viable_count,1)];
f_pareto_end = f_current_gen(pareto_bool_end==1,:);
f_pareto_true_end = [f_pareto_end(:,1), -f_pareto_end(:,2)]; % objective 1 is to be
% minimized while objective 2 is to be maximized
color_end_all = ones(size(f_current_gen,1),1).*color(end); 
color_end = ones(size(f_pareto_true_end,1),1).*color(end);
scatter(f_current_gen(:,1),-f_current_gen(:,2),20,color_end_all)
hold on
scatter(f_pareto_true_end(:,1),f_pareto_true_end(:,2),20,color_end,'filled')
legendarg(end-1,1) = strcat('Generation ', string(n_gen));
legendarg(end,1) = strcat('Generation ', string(n_gen),' - non-dominant');

%current_pareto_bool = plot_pareto_seak(f_current_gen,2);
hold off

xlabel('$\left|\frac{C_{22}}{C_{11}} - c_{target}\right|$','Interpreter','latex','FontSize',15,'FontWeight','bold')
ylabel('$\frac{C_{22}}{v_f}$','Interpreter','latex','FontSize',15,'FontWeight','bold')
legend(legendarg)
title('Pareto Front')

end

