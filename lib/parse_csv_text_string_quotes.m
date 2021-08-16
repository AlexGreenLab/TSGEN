function text_array = parse_csv_text_string_quotes(temp_string)
%function text_array = parse_csv_text_string_quotes(temp_string)

index_list = [];
for c1 = 1:length(temp_string)
    if temp_string(c1) == ','
        index_list(end+1) = c1;
    end
end

index_quote_list = findstr(temp_string,'"');
remove_list = [];
for c1 = 1:2:length(index_quote_list)
    left_index = index_quote_list(c1);
    right_index = index_quote_list(c1+1);
    remove_list = [remove_list,find(left_index <= index_list & right_index >= index_list)];
end
index_list(remove_list) = [];
    

index1 = 1;
for c1 = 1:length(index_list)
    if index1 == index_list(c1)
        text_array{1,c1} = '';
    else
        text_array{1,c1} = temp_string(index1:index_list(c1)-1);
    end
    index1 = index_list(c1)+1;
end
if index1 <= length(temp_string)
    text_array{1,c1+1} = temp_string(index1:end);
end
for c1 = 1:length(text_array)
    temp_string = text_array{c1};
    temp_string(findstr(temp_string,'"')) = [];
    text_array{c1} = temp_string;
end