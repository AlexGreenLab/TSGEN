function output_string = r_revcomp(input_string)
%function output_string = r_revcomp(input_string)

output_string = dna2rna(dna_comp_strand_string(rna2dna(input_string)));
temp = findstr(fliplr(input_string),'N');
output_string(temp) = 'N';