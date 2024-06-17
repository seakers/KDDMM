function index = find_closest_index(array, value)
% Finds the closest index in the array corresponding to the given value

diff_array = array - ones(size(array,1),1).*value;

min_val = min(abs(diff_array));

index = find(abs(diff_array) == min_val,1,'last');

end