function C = SMatrixDecay(Q, times, nbins)
% C = SMatrixDecay(Q, times, nbins)
%
% It computes the decay of correlation between state vectors 
%
% INPUTS: 
% Q: a Qmatrix (rows=times, columns=cells)
% times: the times where to trigger the correlation
% nbins: semiwidth of the histogram
%
% OUTPUTS:
% C: the histogram

% batta 2000
% status: beta

[nQ, c] = size(Q);

times = times(find(times >= nbins+1 & times <= nQ-nbins));
mQ = mean(Q, 1);
Q = Q - repmat(mQ, length(Q), 1);


sQ2 = sqrt(mean(Q.*Q, 2)); % column vector

CT = zeros(2 * nbins+1, length(times));

T0 = Q(times, :);
sQ2_1 = sQ2(times);
warning off
for M = -nbins:nbins
  T1 = Q(times+M,:);
  R = mean(T0 .* T1, 2);
  R = R ./ (sQ2_1 .* sQ2(times+M));
  R(find(isnan(R))) = 0;
  CT(M+nbins+1,:) = R';
end
warning on
C = mean(CT, 2);


















