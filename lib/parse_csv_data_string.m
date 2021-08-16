function data_table = parse_csv_data_string(temp_string)
%function data_table = parse_csv_data_string(temp_string)

index_list = [];
for c1 = 1:length(temp_string)
    if temp_string(c1) == ','
        index_list(end+1) = c1;
    end
end

index1 = 1;
for c1 = 1:length(index_list)
    if index1 == index_list(c1)
        data_table(1,c1) = 0;
    else
        data_table(1,c1) = str2double(temp_string(index1:index_list(c1)-1));
    end
    index1 = index_list(c1)+1;
end
if index1 <= length(temp_string)
    data_table(1,c1+1) = str2double(temp_string(index1:end));
end