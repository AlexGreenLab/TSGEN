

function defect_level = checkUnpairedT(rseq,T)
%function defect_level = checkUnpairedT(rseq,T)	

	dots = char(zeros(1,1000)+'.');
	if length(rseq) == 0
		defect_level = 1;
		disp('Error in checkUnpaired: length(rseq) = 0.');
		return;
	end
	filename = randseq(8);
	fid = fopen([filename,'.in'],'w');
	fprintf(fid,'%s\n%s',rseq,dots(1:length(rseq)));
	fclose(fid);
	a = 1;
	counter = 1;
	while a ~= 0 && counter < 5
		a = system(sprintf('complexdefect -T %d -material rna %s > %s.out',T,filename,filename));
		counter = counter + 1;
	end
	defect_level = loadNupackDefect([filename,'.out'],2);
	if length(defect_level) == 0
		defect_level = 1;
	end
	system(sprintf('rm %s.in %s.out',filename,filename));
end