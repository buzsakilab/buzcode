function M = spkZcorrcoef(S,binSize,is,TT)

%  M = SPKCORRCOEF(S,binSize,is)
%  this function computes correlation coefficients matrix of 
%  spikes in tsa *S*, split in bins of *binSize* for interval *is*

Q = MakeQfromS(S,binSize);
Q = Restrict(Q,is);
M = full(nancorrcoef(Data(Q)));
M(isnan(M))=0;
M = M - diag(diag(M));


%Avoid correlation computation between cells of the same TT

nbCells = length(S);

for i=1:nbCells

	for j=i+1:nbCells
		if (TT(i) == TT(j))
			M(i,j) = 0;
			M(j,i) = 0;
		end
	end

end
