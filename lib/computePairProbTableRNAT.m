function pair_prob_table = computePairProbTableRNAT(rseq,T)
%function pair_prob_table = computePairProbTableRNAT(rseq,T)

	if length(rseq) == 0
		struc_out = char(zeros(1,10)+'(');
		return 
	end
	filename = randseq(8);
	fid = fopen([filename,'.in'],'w');
	fprintf(fid,'%s',rseq);
	fclose(fid);
	a = 1;
	counter = 1;
	while a ~= 0 && counter < 5
		a = system(sprintf('pairs -material rna -T %d %s',T,filename));
		if a ~= 0
			continue;
		end
		pair_prob_table = loadNupackPairProbTable([filename,'.ppairs']);
		if isempty(pair_prob_table)
			fprintf('computePairProbTableRNA: pair_prob_table is empty.\n');
			a = 1;
			continue;
		end
		counter = counter + 1;
	end

	system(sprintf('rm %s.in %s.ppairs',filename,filename));
end
