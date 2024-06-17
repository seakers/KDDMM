function is_diagonal = identify_diagonal_members(sidenum) 
    
    n_total_members = nchoosek(sidenum^2,2);
    is_diagonal = zeros(n_total_members,1);
    member_count = 0;
    
    for i = 1:(sidenum^2)
        closest_top_node = find_closest_top_node(i,sidenum);
        for j = i+1:(sidenum^2)
            if ((j <= closest_top_node) || (rem((j-i),sidenum) == 0))
                is_diagonal(member_count+1,1) = 0;
            else
                is_diagonal(member_count+1,1) = 1;
            end
            member_count = member_count + 1;
        end
    end

end