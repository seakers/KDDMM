% Number of Members Heuristic
% This heuristic function identifies the number of members in each design,
% and scores the design relative to a threshold representative of the
% minimum number of members desired for stable designs
function numMemScore = numMembersHeuristic_2D(CA,threshold,totPosNum)
    numMems = size(CA,1);
    if threshold > (totPosNum/2)
        div = 1/abs(threshold);
    else
        div = 1/abs(threshold - totPosNum);
    end
    numMemScore = round((1-(0.1*div*abs(threshold-numMems))),1);
end