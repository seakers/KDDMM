function top_edge_nodes = get_top_edge_nodes(sidenum)
    reached_right_edge = false;
    top_edge_nodes = zeros(sidenum,1);
    node = sidenum;
    count = 0;
    while(~reached_right_edge)
        if (node > (sidenum^2))
            reached_right_edge = true;
        else
            top_edge_nodes(count+1) = node;
            node = node + sidenum;
            count = count + 1;
        end
    end
end