function fid = writeCSVFileSpectrum(FileName, CSV, Labels)
%function fid = writeCSVFileSpectrum(FileName, CSV, Labels)
fid = fopen(FileName,'w');
for c1 = 1:size(Labels,1)
    for c2 = 1:size(Labels,2)
        if size(Labels{c1,c2},1) == 0
            fprintf(fid, ',');
        else
            fprintf(fid, '%s,', Labels{c1,c2});
        end
    end
    fprintf(fid, '\n');
end

for c1 = 1:size(CSV,1)
    flag = 1;
    for c2 = 1:size(CSV,2)
        if CSV(c1,c2)*0 == 0
            flag = 0;
            break;
        end
    end
    if flag == 1
        continue;
    end
    fprintf(fid,'%12.12f,', CSV(c1,:));
    fprintf(fid, '\n');
end
fclose(fid);
