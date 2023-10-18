function is_repeated_edge_member = identify_repeated_edge_members(sidenum)

    CA_all = get_CA_all(sidenum);
    top_edge_nodes = sidenum:sidenum:sidenum^2;
    member_count = 0;
    is_repeated_edge_member = zeros(size(CA_all,1),1);
    for i = 1:sidenum^2-1
        right_edge = false;
        top_edge = false;
        for j = i+1:sidenum^2
            if i > (sidenum^2-sidenum)
                if j > (sidenum^2-sidenum)
                    is_repeated_edge_member(member_count+1) = 1;
                    right_edge = true;
                end
            end
            if ismember(i,top_edge_nodes)
                if ismember(j,top_edge_nodes)
                    is_repeated_edge_member(member_count+1) = 1;
                    top_edge = true;
                end
            end
            if (~right_edge && ~top_edge)
                is_repeated_edge_member(member_count+1) = 0;
            end
            member_count = member_count + 1;
        end
    end
                    
end