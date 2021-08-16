function line_data = parseNupackEqFileLine(fid)

temp_string = fgetl(fid);
if temp_string(1) == '%'
    line_data = [];
    return;
end
curr_point = 1;
pos_list = findstr(temp_string,char(9));
line_data = [];
for c1 = 1:length(pos_list)
    line_data(1,end+1) = str2num(temp_string(curr_point:pos_list(c1)-1));
    curr_point = pos_list(c1)+1;
end
