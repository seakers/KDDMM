function repeat_member = get_repeated_member(member,sidenum)

    repeat_member = [];
    left_edge_nodes = 1:1:sidenum;
    bottom_edge_nodes = 1:sidenum:(sidenum*sidenum - sidenum)+1;
    top_edge_nodes = sidenum:sidenum:sidenum*sidenum+1;
    right_edge_nodes = (sidenum*sidenum - sidenum)+1:1:sidenum*sidenum;
    
    if (all(ismember(member,left_edge_nodes))) % member is a left edge member
        repeat_member = member + (sidenum*sidenum - sidenum);
    elseif (all(ismember(member,bottom_edge_nodes))) % member is a bottom edge member
        repeat_member = member + sidenum - 1;
    elseif (all(ismember(member,right_edge_nodes))) % member is a right edge member
        repeat_member = member - (sidenum*sidenum - sidenum);
    elseif (all(ismember(member,top_edge_nodes))) % member is a top edge member
        repeat_member = member - (sidenum - 1);
    end
end