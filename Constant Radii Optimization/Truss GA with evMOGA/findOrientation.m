function orientation = findOrientation(p,q,r)
% This function finds the orientation of an ordered triplet (p, q, r)
% The function returns one of three following values:
% 0 --> p, q and r are colinear 
% 1 --> Clockwise 
% 2 --> Counterclockwise 
    val = ((q(2)-p(2))*(r(1)-q(1)))-((q(1)-p(1))*(r(2)-q(2)));
    if val == 0
        orientation = 0;
    elseif val > 0
        orientation = 1;
    else
        orientation = 2;
    end
end