function A = AutoCorrMaze(A)

A = getResource(A,'SpikeData');
nbCells = length(S);

A = getResource(A,'CellNames');
nbCells = length(S);

A = getResource(A,'Sleep1Epoch')
sleep1Epoch = sleep1Epoch{1};

A = getResource(A,'Sleep2Epoch')
sleep2Epoch = sleep2Epoch{1};


resdir = [parent_dir(A), filesep 'XCorr' filesep 'Sleep' filesep 'Auto'];
[p,ds,e] = fileparts(current_dir(A));

for i=1:nbCells

	t = Range(Restrict(S{i},sleep1Epoch));
	[C,B] = CrossCorr(t,t,1,600);
	C(301) = 0;
	
	fh = figure(1),clf
	bar(B,C);
	
	saveas(fh,[resdir filesep ds '_' cellnames{i} 'AutoCorrSleep1_300ms_1msBin'],'png')

end


for i=1:nbCells

	t = Range(Restrict(S{i},sleep2Epoch));
	[C,B] = CrossCorr(t,t,1,600);
	C(301) = 0;
	
	fh = figure(1),clf
	bar(B,C);
	
	saveas(fh,[resdir filesep ds '_' cellnames{i} 'AutoCorrSleep2_300ms_1msBin'],'png')

end