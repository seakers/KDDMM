function complete_bool_array = get_complete_boolean_array(bool_des,sidenum)
    top_nodes = get_top_edge_nodes(sidenum);
    CA_rep = get_repeated_CA(sidenum);
    %n_total_members = nchoosek(sidenum^2,2);
    complete_bool_array = [];
    n_count = 0;
    for i=1:(sidenum^2)
        for j=i+1:(sidenum^2)
            right_edge = false;
            top_edge = false;
            if (i > (sidenum^2 - sidenum))
                if (j > (sidenum^2 - sidenum))
                    repeated_member = [i - ((sidenum-1)*(sidenum)), j - ((sidenum-1)*sidenum)];
                    repeat_index = find_member_index(CA_rep, repeated_member);
                    complete_bool_array = [complete_bool_array;bool_des(repeat_index)];
                    right_edge = true;
                end
            end
            if ismember(i,top_nodes)
                if ismember(j,top_nodes)
                    repeated_member = [i - (sidenum-1), j - (sidenum-1)];
                    repeat_index = find_member_index(CA_rep, repeated_member);
                    complete_bool_array = [complete_bool_array;bool_des(repeat_index)];
                    top_edge = true;
                end
            end
            if (~right_edge && ~top_edge)
                complete_bool_array = [complete_bool_array;bool_des(n_count+1)];
                n_count = n_count + 1;
            end
        end
    end        
end