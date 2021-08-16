function rna_string = dna_matrix2string(rna_mat)
%function dna_string = dna_matrix2string(dna_mat)

rna_string = zeros(1,length(rna_mat));
for c1 = 1:length(rna_mat)
    switch rna_mat(c1)
        case 0
            rna_string(c1) = 'A';
            continue;
        case 3
            rna_string(c1) = 'T';
            continue;
        case 1
            rna_string(c1) = 'G';
            continue;
        case 2
            rna_string(c1) = 'C';
            continue;
    end
end
rna_string = char(rna_string);