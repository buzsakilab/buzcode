function [outR, outID] = makeBayesWeightedCorrBatch1(Pr, bID, varargin)
%function [outR, outID] = makeBayesWeightedCorrBatch1(Pr, bID, preComputedQ (optional))
% inputs 
% Pr: [nTimeBin X nSpatialBin] matrix of posterior probabilities
%       (each row should sum to 1)
% bID = 'bin Id's', [nTimeBin X 1] vector
%       if only 1 event is inputed then bID is a vector ones, 
%       the weighted correlations will be calculated on each unique value
%       in bID
% Q (optional input) - Q is a linearized list of space X time-bin indeces
%       can be precomputed to save on computational time.

outID = unique(bID);
h = histc(bID, outID);
if isempty(varargin)
    Q = makeQForWeightedCorr(h, size(Pr, 2));
else
    Q = varargin{1};
end

outR = zeros(length(outID), 1);
uH = unique(h);
nB = size(Pr, 2);
for i = 1:length(uH)
    k1 = ismember(bID, outID(h == uH(i)));
    g3 = reshape(permute(reshape(Pr(k1 == 1, :), [uH(i), sum(h == uH(i)), nB]), [1, 3, 2]), [uH(i)*nB, sum(h == uH(i))]);
   
    outR(h == uH(i)) = makeWeightedCorrBatchIn(Q{uH(i)}, g3);
end

end

function out = makeWeightedCorrBatchIn(xy, w);
%function out = makeWeightedCorr1(xy, w);

x = repmat(xy(:, 1), 1, size(w, 2));
y = repmat(xy(:, 2), 1,  size(w, 2));

mx = repmat(sum(w.*x, 1)./sum(w, 1), size(w, 1), 1);
my = repmat(sum(w.*y, 1)./sum(w, 1), size(w, 1), 1);

covxy = sum(w.*(x - mx).*(y - my))./sum(w);
covxx = sum(w.*((x - mx).^2))./sum(w);
covyy = sum(w.*((y - my).^2))./sum(w);

out = (covxy./(sqrt(covyy.*covxx)))';
end
