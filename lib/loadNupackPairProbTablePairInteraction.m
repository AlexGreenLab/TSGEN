function [pair_prob_table_set,combo_set] = loadNupackPairProbTablePairInteraction(FileName)

combo_set = [1,0;0,1;2,0;1,1;0,2];
out_seq = 1;
fid = fopen(FileName,'r');
if fid == -1
    return;
end
temp_string = fgetl(fid);
out_seq = 1;
while ischar(temp_string) && isempty(temp_string) ~= 1
    if temp_string(1) ~= '%'
        break;
    end
    temp_string = fgetl(fid);
end

for c0 = 1:5
    temp_string = fgetl(fid);
    temp_string = fgetl(fid);
    seq_length = str2num(fgetl(fid));
    pair_prob_table = zeros(seq_length,seq_length);
    temp_string = fgetl(fid);
    while ischar(temp_string) && isempty(temp_string) ~= 1
        if temp_string(1) == '%'
            break;
        end
        tab_pos = findstr(temp_string,char(9));
        if numel(tab_pos) == 1
            pair_prob_table = [];
            fclose(fid);
            return;
        end
        index1 = str2num(temp_string(1:tab_pos(1)-1));
        index2 = str2num(temp_string(tab_pos(1)+1:tab_pos(2)-1));
        if index2 > seq_length
            index2 = index1;
        end
        if index2 < index1
            disp('index2 < index1');
        end
        curr_prob = str2num(temp_string(tab_pos(2)+1:end));
        pair_prob_table(index1,index2) = curr_prob;
        pair_prob_table(index2,index1) = curr_prob;
        temp_string = fgetl(fid);
    end
    temp_string = fgetl(fid);
    for c1 = 1:size(pair_prob_table,1)
        if pair_prob_table(c1,c1) == 0 && sum(pair_prob_table(c1,:)) < 1
            pair_prob_table(c1,c1) = 1 - sum(pair_prob_table(c1,:));
        end
    end
    pair_prob_table_set{c0,1} = pair_prob_table;
end

fclose(fid);
