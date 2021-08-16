function defect_level = checkDefectArbT(rseq,rstruc,T)
%function defect_level = checkDefectArbT(rseq,rstruc,T)

    filename = randseq(16);
	fid = fopen([filename,'.in'],'w');
	fprintf(fid,'%s\n%s',rseq,rstruc);
	fclose(fid);
	a = 1;
	counter = 1;

	while a ~= 0 && counter < 50
		a = system(sprintf('complexdefect -T %d -material rna %s > %s.out',T,filename,filename));
		counter = counter + 1;
	end

	defect_level = loadNupackDefect([filename,'.out'],2);
	system(sprintf('rm %s.in %s.out',filename,filename));

end