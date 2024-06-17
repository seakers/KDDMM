function member_index = find_member_index(CA_full, member)
    member_index = 0;
    for i=1:size(CA_full,1)
        if (CA_full(i,1) == member(1))
            if (CA_full(i,2) == member(2))
                member_index = i;
            end
        end
    end
end