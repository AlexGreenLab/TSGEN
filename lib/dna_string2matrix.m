function rna_mat = dna_string2matrix(rna_string)
%function dna_mat = dna_string2matrix(dna_string)

rna_mat = zeros(1,length(rna_string));
for c1 = 1:length(rna_string)
    switch rna_string(c1)
        case 'a'
            rna_mat(c1) = 0;
            continue;
        case 't'
            rna_mat(c1) = 3;
            continue;
        case 'g'
            rna_mat(c1) = 1;
            continue;
        case 'c'
            rna_mat(c1) = 2;
            continue;
        case 'A'
            rna_mat(c1) = 0;
            continue;
        case 'T'
            rna_mat(c1) = 3;
            continue;
        case 'G'
            rna_mat(c1) = 1;
            continue;
        case 'C'
            rna_mat(c1) = 2;
            continue;
    end
end