function contains = map_contains_design(design_map, design) 
    contains = false;
    for j = keys(design_map)
        key = j{1};
        if isequal(design_map(key),design)
            contains = true;
            break
        end
    end
end