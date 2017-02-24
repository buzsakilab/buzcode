function M = spkcorrcoef(Q,TT)

%  Y = spkcorrcoef(S,binSize,is,TT)
%  this function computes correlation coefficients matrix of 
%  spikes in tsa *S* whith Names TT, split in bins of *binSize* for interval *is*

M = nancorrcoef(full(Data(Q))); %useless to use zscore here
M(isnan(M))=0;
M = M - diag(diag(M));

%Avoid correlation computation between cells of the same TT

nbCells = size(M,1);

for i=1:nbCells

	for j=i+1:nbCells
		if (TT(i) == TT(j))
			M(i,j) = 0;
			M(j,i) = 0;
		end
	end

end
