function M = spkcorrcoef(S,binSize,is)

%  Y = spkcorrcoef(S,binSize,is,TT)
%  this function computes correlation coefficients matrix of 
%  spikes in tsa *S* whith Names TT, split in bins of *binSize* for interval *is*

Q = MakeQfromS(S,binSize);
Q = Restrict(Q,is);
M = full(nancorrcoef(Data(Q))); %useless to use zscore here
M(isnan(M))=0;
M = M - diag(diag(M));

