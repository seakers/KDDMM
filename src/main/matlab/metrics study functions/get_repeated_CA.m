function CA_rep = get_repeated_CA(sidenum)
    top_nodes = get_top_edge_nodes(sidenum);
    n_member = 0;
    n_total_members = nchoosek(sidenum^2,2);
    CA_rep_all = zeros(n_total_members,2);
    for i = 1:(sidenum^2)
        for j = i+1:(sidenum^2)
            % Check (and don't include) if member is on the right edge
            if (i > (sidenum^2 - sidenum))
                if (j > (sidenum^2 - sidenum))
                    continue
                end
            end
            % Check (and don't include) if member is on the top edge
            if ismember(i,top_nodes)
                if ismember(j,top_nodes)
                    continue
                end
            end
            CA_rep_all(n_member+1,:) = [i,j];
            n_member = n_member + 1;
        end
    end
    CA_rep = CA_rep_all(1:n_member,:);
end