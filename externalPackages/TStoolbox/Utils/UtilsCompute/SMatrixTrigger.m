function S = SmatrixTrigger(Q, times, nbins)
% C = SMatrixDecay(Q, times, nbins)
%
% It computes the average of the S matrix triggered by the events times
%
% INPUTS: 
% Q: a Qmatrix (rows=times, columns=cells)
% times: the times where to trigger the Smatrix
% nbins: semiwidth of the matrix
%
% OUTPUTS:
% C: the histogram

% batta 2000
% status: beta

[nQ, c] = size(Q);

times = times(find(times >= nbins+1 & times <= nQ-nbins));
mQ = mean(Q,  1);
Q = Q - repmat(mQ, length(Q), 1);

S = zeros(2 * nbins + 1, 2 * nbins + 1);

for i = times
  S = S + corrcoef(Q(i-nbins:i+nbins, :)');
end

S  = S / length(times);
