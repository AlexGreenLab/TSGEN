function dna_comp_mat = dna_comp_strand(dna_mat)
%function dna_comp_mat = dna_comp_strand(dna_mat)

dna_comp_mat = 3 - dna_mat;
dna_comp_mat = fliplr(dna_comp_mat);