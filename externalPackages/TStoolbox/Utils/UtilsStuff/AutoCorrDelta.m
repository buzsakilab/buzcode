function A = AutoCorrMaze(A)

A = getResource(A,'SpikeData');
nbCells = length(S);

A = getResource(A,'CellNames');
nbCells = length(S);

A = getResource(A,'Sleep1DeltaEpochs')
sleep1DeltaEpochs = sleep1DeltaEpochs{1};

A = getResource(A,'Sleep2DeltaEpochs')
sleep2DeltaEpochs = sleep2DeltaEpochs{1};



resdir = [parent_dir(A), filesep 'XCorr' filesep 'Sleep' filesep 'Auto'];
[p,ds,e] = fileparts(current_dir(A));

for i=1:nbCells

	t = Range(Restrict(S{i},sleep1DeltaEpochs));
	[C,B] = CrossCorr(t,t,1,2000);
	C(1001) = 0;
	
	fh = figure(1),clf
	bar(B,C);
	
	saveas(fh,[resdir filesep ds '_' cellnames{i} 'AutoCorrDeltaSleep1_1s_1msBin'],'png')

end


for i=1:nbCells

	t = Range(Restrict(S{i},sleep2DeltaEpochs));
	[C,B] = CrossCorr(t,t,1,2000);
	C(1001) = 0;
	
	fh = figure(1),clf
	bar(B,C);
	
	saveas(fh,[resdir filesep ds '_' cellnames{i} 'AutoCorrDeltaSleep2_1s_1msBin'],'png')

end