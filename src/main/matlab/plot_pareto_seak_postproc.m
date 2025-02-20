function [] = plot_pareto_seak_postproc(f_vals, NC, CA_all, design_array, feas_array, stab_array, plot_case, n_objs, marker)
% This function uses paretofront.m, paretofront.c and paretofront.mexw64 to 
% determine the designs in the pareto front and plots them. (add figure
% statement in main script)
% Inputs: Objective values for all evaluated designs (n x 2)
%         Number of objectives n_objs

pareto_bool = paretofront(f_vals);

f_pareto = f_vals(pareto_bool==1,:);

if (n_objs == 2)
    %%% For 2 objectives
    f_pareto_true = [f_pareto(:,1), -f_pareto(:,2)]; % objective 1 is to be
    % minimized while objective 2 is to be maximized  

    plot(f_pareto_true(:,1),f_pareto_true(:,2),marker,'ButtonDownFcn',{@plot_pareto_callback,NC,CA_all,design_array,feas_array,stab_array,f_vals})
    switch plot_case
        case 'true_objectives'
            xlabel('$\left|\frac{C_{22}}{C_{11}} - c_{target}\right|$','Interpreter','latex','FontSize',15,'FontWeight','bold')
            ylabel('$\frac{C_{22}}{v_f}$','Interpreter','latex','FontSize',15,'FontWeight','bold')
        case 'pen_objectives'
            xlabel('$\frac{\left|\frac{C_{22}}{C_{11}} - c_{target}\right|}{15} - k_{pen}\frac{log(g_{feas}) + log(g_{stab})}{2}$','Interpreter','latex','FontSize',12,'FontWeight','bold')
            ylabel('$\frac{\frac{C_{22}}{v_f}}{8500} + k_{pen}\frac{log(g_{feas}) + log(g_{stab})}{2}$','Interpreter','latex','FontSize',12,'FontWeight','bold')
    end
    title('Pareto Front')

elseif (n_objs == 3)
    %%% For 3 objectives
    f_pareto_true = [f_pareto(:,1), -f_pareto(:,2), -f_pareto(:,3)]; % objective 1 is to be
    % minimized while objectives 2 and 3 are to be maximized  

    plot3(f_pareto_true(:,1),f_pareto_true(:,2),f_pareto_true(:,3),'*')
    xlabel('$\frac{C_{22}}{C_{11}} - c_{target}$','Interpreter','latex','FontSize',15,'FontWeight','bold')
    ylabel('$\frac{C_{22}}{v_f}$','Interpreter','latex','FontSize',15,'FontWeight','bold')
    zlabel('$g_{feas}+g_{stab}$','Interpreter','latex','FontSize',15,'FontWeight','bold')
    title('Pareto Front')

end

end

