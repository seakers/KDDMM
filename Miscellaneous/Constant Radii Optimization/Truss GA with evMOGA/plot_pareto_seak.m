function [] = plot_pareto_seak(f_vals, n_objs)
% This function uses paretofront.m, paretofront.c and paretofront.mexw64 to 
% determine the designs in the pareto front and plots them.
% Inputs: Objective values for all evaluated designs (n x 2)
%         Number of objectives n_objs

pareto_bool = paretofront(f_vals);

f_pareto = f_vals(pareto_bool==1,:);

if (n_objs == 2)
    %%% For 2 objectives
    f_pareto_true = [f_pareto(:,1), -f_pareto(:,2)]; % objective 1 is to be
    % minimized while objective 2 is to be maximized  

    figure
    plot(f_pareto_true(:,1),f_pareto_true(:,2),'*')
    xlabel('$\left|\frac{C_{22}}{C_{11}} - c_{target}\right|$','Interpreter','latex')
    ylabel('$\frac{C_{22}}{v_f}$','Interpreter','latex')
    title('Pareto Front')

elseif (n_objs == 3)
    %%% For 3 objectives
    f_pareto_true = [f_pareto(:,1), -f_pareto(:,2), -f_pareto(:,3)]; % objective 1 is to be
    % minimized while objectives 2 and 3 are to be maximized  

    figure
    plot3(f_pareto_true(:,1),f_pareto_true(:,2),f_pareto_true(:,3),'*')
    xlabel('$\frac{C_{22}}{C_{11}} - c_{target}$','Interpreter','latex')
    ylabel('$\frac{C_{22}}{v_f}$','Interpreter','latex')
    zlabel('$g_{feas}+g_{stab}$','Interpreter','latex')
    title('Pareto Front')

end

end

