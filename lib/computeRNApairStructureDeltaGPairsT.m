function [conc1,conc2,conc_dimer,pair_structure,deltaG,pair_prob_table_set] = computeRNApairStructureDeltaGPairsT(seq1,seq2,T)
% function [conc1,conc2,conc_dimer,pair_structure,deltaG,pair_prob_table_set] = computeRNApairStructureDeltaGPairsT(seq1,seq2,T)

	if length(seq1) == 0
		defect_level = 1;
		disp('Error in computeDNAprimerpairs: length(seq1) = 0.');
		return;
	end
	if length(seq2) == 0
		defect_level = 1;
		disp('Error in computeDNAprimerpairs: length(seq2) = 0.');
		return;
	end

	filename = randseq(8);
	fid = fopen([filename,'.in'],'w');
	fprintf(fid,'2\n%s\n%s\n2',seq1,seq2);
	fclose(fid);
	fid = fopen([filename,'.con'],'w');
	fprintf(fid,'1e-7\n1e-7');
	fclose(fid);
	a = 1;
	counter = 1;
	while a ~= 0 && counter < 5
		[a,b] = system(sprintf('complexes -material rna -T %d -ordered -pairs -mfe %s',T,filename));
		[a,b] = system(sprintf('concentrations -T %d -ordered -pairs -mfe %s',T,filename));
		counter = counter + 1;
	end
	conc1 = loadNupackEqFileLineOrder([filename,'.eq'],1,1)./1e-7;
	conc2 = loadNupackEqFileLineOrder([filename,'.eq'],2,1)./1e-7;
	conc_dimer = loadNupackEqFileLineOrder([filename,'.eq'],4,1)./1e-7;

	conc1 = conc1(end);
	conc2 = conc2(end);
	conc_dimer = conc_dimer(end);
	pair_structure{1} = loadNupackMFEFileStructureOrder([filename,'.ocx-mfe'], 1,1);
	pair_structure{2} = loadNupackMFEFileStructureOrder([filename,'.ocx-mfe'], 2,1);
	pair_structure{3} = loadNupackMFEFileStructureOrder([filename,'.ocx-mfe'], 4,1);
	[discard,deltaG(1)] = loadNupackMFEdeltaG([filename,'.ocx-mfe'], 1,1);
	[discard,deltaG(2)] = loadNupackMFEdeltaG([filename,'.ocx-mfe'], 2,1);
	[discard,deltaG(3)] = loadNupackMFEdeltaG([filename,'.ocx-mfe'], 4,1);
	[pair_prob_table_set,~] = loadNupackPairProbTablePairInteraction([filename,'.ocx-ppairs']);
	[a,b] = system(sprintf('rm %s.*',filename));

end

