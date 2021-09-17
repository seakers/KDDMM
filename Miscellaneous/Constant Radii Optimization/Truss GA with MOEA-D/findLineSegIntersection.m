function intersect = findLineSegIntersection(p1,q1,p2,q2)
% This boolean function determines whether two line segments intersect,
% given their endpoints as inputs
    if (findOrientation(p1,q1,p2) ~= findOrientation(p1,q1,q2))&&...
            (findOrientation(p2,q2,p1) ~= findOrientation(p2,q2,q1)) 
        if isequal(p1,p2) || isequal(q1,q2) || ...
                isequal(p1,q2) || isequal(q1,p2)
            intersect = false; 
        else
            intersect = true;
        end
    else
        intersect = false;
    end
end