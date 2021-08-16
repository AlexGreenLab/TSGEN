function array = loadCSVabsFileRawtoStructureQuotes(FileName, PathName)
%function array = loadCSVabsFileRawtoStructureQuotes(FileName, PathName)

TotalName = [PathName,FileName];
FileID = fopen(TotalName);

array = [];
counter = 1;
temp_string = fgetl(FileID);
while ~feof(FileID) || length(temp_string) ~= 0
    if isnumeric(temp_string) == 1
        break;
    end
    if ~isempty(temp_string)
        line1 = parse_csv_text_string_quotes(temp_string);
    else
        line1 = '';
    end
    for c1 = 1:length(line1)
        array{counter,c1} = line1{c1};
    end
    temp_string = fgetl(FileID);
    counter = counter + 1;
end
fclose(FileID);
if size(array{1,1},2) >= 3 && issame(array{1,1}(1:3),'﻿')
    array{1,1} = array{1,1}(4:end);
end
