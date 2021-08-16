function dna = rna2dna_B(rna)
%function dna = rna2dna_B(rna)
% Converts any uracil nucleotides in an RNA sequence into thymine (U-->T).
%
% DNA is returned in the same format as RNA, so if RNA is an integer
% sequence then so is DNA.
% Does not warn user if DNA contains U or other bases besides 'ATCG'

% If the input is a structure then extract the Sequence data.
if isstruct(rna)
    rna = seqfromstruct(rna);
end
if ~ischar(rna)
    dna = rna;
    return
end

dna = strrep(rna,'U','T');
dna = strrep(dna,'u','t');




