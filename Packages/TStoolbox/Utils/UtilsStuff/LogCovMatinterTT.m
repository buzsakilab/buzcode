function A = LogCovMatinterTT(A)

%  z-transform of the log of the column vector X. If X is a matrix, z-transform of each vector is computed


A = getResource(A,'CellNames',datasets{day});
A = getResource(A,'MazeEpoch',datasets{day});

A = registerResource(A, 'SpwTrigCorrMS1_PCA', 'tsdArray', {[],1}, ...
    'spwTrigCorrMS1_PCA', ...
    'spwTrigCorrMS sleep1/maze correlation of PC');


Q = MakeQfromS(S,binSize);
Q = Restrict(Q,is);
M = full(nancorrcoef(zscore(log(Data(Q)+10^-6))));
M(isnan(M))=0;
M = M - diag(diag(M));

mazeEpoch = mazeEpoch{1};
binSpk = MakeQfromS(S,binSizeSleep)
cM = spkZcorrcoef(S,binSizeMaze,mazeEpoch);

nbCells = size(zFiringS2,2);

for i=1:nbCells

	for j=i+1:nbCells
		cell1 = cellnames{i};
		cell2 = cellnames{j};
		if (cell1(3) == cell2(3))
			cM(i,j) = 0;
			cM(j,i) = 0;
		end
	end

end


