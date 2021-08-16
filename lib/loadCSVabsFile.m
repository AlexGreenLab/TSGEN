function [CSV,Labels] = loadCSVabsFile(FileName, PathName)
%function [CSV,Labels] = loadCSVabsFile(FileName, PathName)

TotalName = [PathName,FileName];
FileID = fopen(TotalName);

Labels = [];
temp_string = fgetl(FileID);
line1 = parse_csv_text_string(temp_string);
for c1 = 1:size(line1,2)
    if c1 <= size(line1,2)
        Labels{1,c1} = line1{1,c1};
    else
        Labels{1,c1} = '';
    end
end

CSV = [];
temp_string = fgetl(FileID);
while length(temp_string) ~= 0
    data_line = parse_csv_data_string(temp_string);
    if length(data_line) == 0
        break;
    end
    CSV(end+1,:) = data_line;
    temp_string = fgetl(FileID);
end
fclose(FileID);
