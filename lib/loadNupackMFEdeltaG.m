function [line_data,deltaG] = loadNupackMFEdeltaG(file_name, complex_num, order)
%function [line_data,deltaG] = loadNupackMFEdeltaG(file_name, complex_num, order)


fid = fopen(file_name,'r');
%move to complex data
search_string = [' complex',num2str(complex_num),'-order',num2str(order)];
temp_string = fgetl(fid);
while ischar(temp_string)
    if length(findstr(temp_string,search_string)) ~= 0
        break;
    end
    temp_string = fgetl(fid);
end
if isnumeric(temp_string)
    fclose(fid);
    fid = fopen(file_name,'r');
    search_string = [' composition',num2str(complex_num),'-ordering',num2str(order)];
    temp_string = fgetl(fid);
    while ischar(temp_string)
        if length(findstr(temp_string,search_string)) ~= 0
            break;
        end
        temp_string = fgetl(fid);
    end
    temp_string = fgetl(fid);
    deltaG = str2num(fgetl(fid));
    line_data = fgetl(fid);
    fclose(fid);
else
    temp_string = fgetl(fid);
    deltaG = str2num(fgetl(fid));
    line_data = fgetl(fid);
    fclose(fid);
end
