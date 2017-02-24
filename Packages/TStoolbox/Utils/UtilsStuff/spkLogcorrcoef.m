function M = spkLogcorrcoef(S,binSize,is,cellnames)

%  Y = spkLogTransform(X)
%  z-transform of the log of the column vector X. If X is a matrix, z-transform of each vector is computed

Q = MakeQfromS(S,binSize);
Q = Restrict(Q,is);
M = full(nancorrcoef(zscore(log(Data(Q)+10^-6))));
M(isnan(M))=0;
M = M - diag(diag(M));

%Avoid correlation computation between cells of the same TT

nbCells = length(S);

for i=1:nbCells

	for j=i+1:nbCells
		cell1 = cellnames{i};
		cell2 = cellnames{j};
		if (cell1(3) == cell2(3))
			M(i,j) = 0;
			M(j,i) = 0;
		end
	end

end
