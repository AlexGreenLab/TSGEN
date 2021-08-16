function line_data = loadNupackEqFileLine(file_name, index1,index2)
%function line_data = loadNupackEqFileLine(file_name, index1,index2)

fid = fopen(file_name,'r');

%move to the first line of output data
temp_string = fgetl(fid);
while ischar(temp_string)
    if length(findstr(temp_string,' T =')) ~= 0
        break;
    end
    temp_string = fgetl(fid);
end
counter = 1;
line_data = parseNupackEqFileLine(fid);
while (length(line_data) == 0 && counter < 50000)
    line_data = parseNupackEqFileLine(fid);
    counter = counter + 1;
end
while length(line_data) ~= 0
    if line_data(1) == index1 && line_data(2) == index2
        break;
    end
    line_data = parseNupackEqFileLine(fid);
end
fclose(fid);
