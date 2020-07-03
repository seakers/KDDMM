function [feas] = truss_feasibility(tr1, tr2, NC)
% Check whether the truss choices are feasible
% Inputs: Endpoints for truss 1 tr1 = [p1_1, p2_1]
%         Endpoints for truss 2 tr2 = [p1_2, p2_2]
%         Nodal Connectivity Matrix NC (9 x 2)
% Output: feasibility boolean feas

%x1_1 = NC(tr1(1),1);
y1_1 = NC(tr1(1),2);
x2_1 = NC(tr1(2),1);
y2_1 = NC(tr1(2),2);

x1_2 = NC(tr2(1),1);
y1_2 = NC(tr2(1),2);
%x2_2 = NC(tr2(2),1);
y2_2 = NC(tr2(2),2);
    
feas = false;
if (x1_2 >= x2_1)
    feas = true;
elseif (x1_2 < x2_1)
    if (min([y1_2,y2_2]) >= max([y1_1,y2_1]))
        feas = true;
    end
end
end
