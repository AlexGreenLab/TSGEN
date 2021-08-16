function [deltaG,rna_struc] = computeRNAdeltaGT(seq,T)
	%function [deltaG,rna_struc] = computeRNAdeltaGT(seq,T)

	if length(seq) == 0
		defect_level = 1;
		disp('Error in checkUnpaired: length(seq) = 0.');
		return;
	end
	filename = randseq(8);
	fid = fopen([filename,'.in'],'w');
	fprintf(fid,'1\n%s\n1',seq);
	fclose(fid);
	a = 1;
	counter = 1;
	while a ~= 0 && counter < 5
		[a,b] = system(sprintf('mfe -T %d -multi -material rna %s',T,filename));
		counter = counter + 1;
	end
	[deltaG,rna_struc] = loadNupackMFEsingleStrand([filename,'.mfe']);
	system(sprintf('rm %s.*',filename));
end
