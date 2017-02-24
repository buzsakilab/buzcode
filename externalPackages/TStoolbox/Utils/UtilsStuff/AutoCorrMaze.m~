function A = AutoCorrMaze(A)

A = getResource(A,'MazeEpoch');
mazeEpoch = mazeEpoch{1};

A = getResource(A,'SpikeData');
nbCells = length(S);

A = getResource(A,'CellNames');

resdir = [parent_dir(A), filesep 'XCorr' filesep 'Maze' filesep 'Auto'];
[p,ds,e] = fileparts(current_dir(A));

for i=1:nbCells

	t = Range(Restrict(S{i},mazeEpoch));
	[C,B] = CrossCorr(t,t,1,600);
	C(301) = 0;
	
	fh = figure(1),clf
	bar(B,C);
	
	saveas(fh,[resdir filesep ds '_' cellnames{i} 'AutoCorrMaze300ms_1msBin'],'png')

end
