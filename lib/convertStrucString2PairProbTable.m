% convertStrucString2PairProbTable returns the result: 
% pair_probe_table, whose diagonal shows probability of unpaired bases. Input
% should use N to specify bases where pairing is not known or considered, 
% string must still have equal number of "(" and ")"

function pair_prob_table = convertStrucString2PairProbTable(struc_string)
%function pair_prob_table = convertStrucString2PairProbTable(struc_string)

%use stack to store indices
seq_length = length(struc_string);
pair_prob_table = zeros(seq_length,seq_length);

index_stack = [];
for c1 = 1:seq_length
    if struc_string(c1) == '.'
        pair_prob_table(c1,c1) = 1;
    elseif struc_string(c1) == '(';
        index_stack(end+1,1) = c1;
    elseif struc_string(c1) == ')';
        if isempty(index_stack)
            disp('Error: struc_string has uneven number of brackets.');
            return;
        end
        index1 = index_stack(end);
        index_stack(end) = [];
        pair_prob_table(index1,c1) = 1;
        pair_prob_table(c1,index1) = 1;
    end
end
        
        
