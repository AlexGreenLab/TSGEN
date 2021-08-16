function rna = dna2rna_B(dna)
%function rna = dna2rna(dna)
% Converts any thymine nucleotides in a DNA sequence into uracil (T-->U).
%
% RNA is returned in the same format as DNA, so if DNA is an integer
% sequence then so is RNA.
% Does not warn user if DNA contains U or other bases besides 'ATCG'


% If the input is a structure then extract the Sequence data.
if isstruct(dna)
    dna = seqfromstruct(dna);
end

if ~ischar(dna)
    rna = dna;
    return
end


rna = strrep(dna,'T','U');
rna = strrep(rna,'t','u');

