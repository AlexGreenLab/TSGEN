function dna_comp_string = d_revcomp(dna_string)
%function dna_comp_string = d_revcomp(dna_string)

if (length(findstr(dna_string,'u')) + length(findstr(dna_string,'U'))) > 0
    dna_string = rna2dna(dna_string);
end
dna_mat = dna_string2matrix(dna_string);
dna_comp_mat = dna_comp_strand(dna_mat);
dna_comp_string = dna_matrix2string(dna_comp_mat);
temp = findstr(fliplr(dna_string),'N');
dna_comp_string(temp) = 'N';