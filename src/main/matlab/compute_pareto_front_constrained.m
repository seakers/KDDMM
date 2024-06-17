% Determine indices of pareto front designs from the given population given constraints
% NOTE: Assumes the minimization of all objectives

function ispareto = compute_pareto_front_constrained(pop_objs, pop_aggr_constrs)
	% pop_objs = array of size n_des x n_objs, pop_aggr_constrs = array of size n_des (sum of absolute value of constraints)
	pop_size = size(pop_objs, 1);
	n_objs = size(pop_objs, 2);
	
	dom_counter = zeros(pop_size, 1);
	
	for i = 1:pop_size
		for j = i+1:pop_size
			% Compare constraint satisfaction
			if pop_aggr_constrs(i) > pop_aggr_constrs(j)
				dom_counter(i) = dom_counter(i) + 1;
			elseif pop_aggr_constrs(j) > pop_aggr_constrs(i)
				dom_counter(j) = dom_counter(j) + 1;
			else
				% In case of ties, compare objective minimization
				dominate = zeros(n_objs, 1);
				for k = 1:n_objs
					if pop_objs(i,k) > pop_objs(j,k)
						dominate(k) = 1;
					elseif pop_objs(i,k) < pop_objs(j,k)
						dominate(k) = -1;
					end
				end
				if (ismember(1, dominate) && ~ismember(-1, dominate))
					dom_counter(i) = dom_counter(i) + 1;
				elseif (ismember(-1, dominate) && ~ismember(1, dominate))
					dom_counter(j) = dom_counter(j) + 1;
				end
			end
		end
	end
	
	ispareto = dom_counter == 0;
    
    %%% Testing
    % pop_objs = [1,7; 3,8; 3,6; 4,6; 2,6; 3,5; 5,5; 4,4; 5,3; 4,8];
    % pop_constr_aggr = [0.5; 0.5; 0; 0; 0.75; 0.6; 0.3; 0.25; 0; 0.8];
    % ispareto = compute_pareto_front_constrained(pop_objs, pop_constr_aggr);
    % pareto_objs = pop_objs(ispareto, :);
    %%% Expeected answer : [3,6; 5,3]
	
end
					
				
				