function x_des = get_binary_array_from_bitstring(des_string)
	x_des = zeros(strlength(des_string),1);
	for i = 1:strlength(des_string)
		x_des(i,1) = str2double(des_string(i));
	end
end