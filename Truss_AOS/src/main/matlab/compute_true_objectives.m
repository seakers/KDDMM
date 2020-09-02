function [f_true] = compute_true_objectives(csv_data, pop_size, fibre_stiffness)
% This function computes the true objectives given the penalized objectives
f_penalized = csv_data{:,1:2};
feas_scores = csv_data{:,3}; 
stab_scores = csv_data{:,4};
f_true = zeros(pop_size, 2);
pen_fac = 1;
if fibre_stiffness
    pen_fac = 1.5;
end
for i = 1:pop_size
    penalty = (log10(abs(feas_scores(i))) + log10(abs(stab_scores(i))))/2;
    f_true(i,:) = [15*(f_penalized(i,1) + pen_fac*penalty), 85000*(f_penalized(i,2) + pen_fac*penalty)];
end

end

