function M = spkcorrcoefAll(Q,TT)

%  Y = spkcorrcoef(S,binSize,is,TT)
%  this function computes correlation coefficients matrix of 
%  spikes in tsa *S* whith Names TT, split in bins of *binSize* for interval *is*

M = nancorrcoef(full(Data(Q))); %useless to use zscore here
M(isnan(M))=0;
M = M - diag(diag(M));
