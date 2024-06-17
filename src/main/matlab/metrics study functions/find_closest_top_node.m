function closest_top_node = find_closest_top_node(node,sidenum)
    
    top_reached = false;
    current_node = node;
    while ~top_reached
        if (rem(current_node,sidenum) == 0)
            top_reached = true;
        else
            current_node = current_node + 1;
        end
    end

    closest_top_node = current_node;
end