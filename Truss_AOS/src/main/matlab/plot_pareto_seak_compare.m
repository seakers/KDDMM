function [] = plot_pareto_seak_compare(f_vals_eps,f_vals_aos)
% This function plots the pareto fronts using epsilon MOEA and  AOS in the
% same plot for comparison 

pareto_bool_eps = paretofront(f_vals_eps);
pareto_bool_aos = paretofront(f_vals_aos);

f_pareto_eps = f_vals_eps(pareto_bool_eps==1,:);
f_pareto_aos = f_vals_aos(pareto_bool_aos==1,:);

f_pareto_true_eps = [f_pareto_eps(:,1), -f_pareto_eps(:,2)]; 
f_pareto_true_aos = [f_pareto_aos(:,1), -f_pareto_aos(:,2)]; % objective 1 is to be
% minimized while objective 2 is to be maximized  

figure
plot(f_pareto_true_eps(:,1),f_pareto_true_eps(:,2),'*b')
hold on
plot(f_pareto_true_aos(:,1),f_pareto_true_aos(:,2),'*r')
hold off
xlabel('$\left|\frac{C_{22}}{C_{11}} - c_{target}\right|$','Interpreter','latex','FontSize',15,'FontWeight','bold')
ylabel('$\frac{C_{22}}{v_f}$','Interpreter','latex','FontSize',15,'FontWeight','bold')
legend('e-MOEA','AOS','Location','Best')
title('Pareto Front Comparison b/w e-MOEA and AOS')

end

