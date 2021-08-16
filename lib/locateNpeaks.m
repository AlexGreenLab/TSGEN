function peak_table = locateNpeaks(N,pos,height,discard_window)
%function peak_table = locateNpeaks(N,pos,height,discard_window)

[pos,indices] = sort(pos);
height = height(indices);
peak_table = [];
while ~isempty(pos) && N > 0
    [max_height,max_index] = max(height);
    max_pos = pos(max_index);
    remove_list = [max_index];
    peak_table(end+1,:) = [max_pos,max_height];
    for c1 = max_index:length(pos)
        if pos(c1) >= max_pos-discard_window && pos(c1) <= max_pos+discard_window
            remove_list(end+1,1) = c1;
        elseif pos(c1) > max_pos+discard_window
            break;
        end        
    end
    for c1 = max_index:-1:1
        if pos(c1) >= max_pos-discard_window && pos(c1) <= max_pos+discard_window
            remove_list(end+1,1) = c1;
        elseif pos(c1) < max_pos-discard_window
            break;
        end        
    end
    pos(remove_list) = [];
    height(remove_list) = [];
    N = N - 1;
end