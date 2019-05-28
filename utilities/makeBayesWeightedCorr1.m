function [outR, outID] = makeBayesWeightedCorr1(Pr, bID, varargin)
%function [outR, outID] = makeBayesWeightedCorr1(Pr, bID, preComputedQ (optional))

outID = unique(bID);
h = histc(bID, outID);
if isempty(varargin)
    Q = makeQForWeightedCorr(h, size(Pr, 2));
else
    Q = varargin{1};
end

outR = zeros(length(outID), 1);
for i = 1:length(outR)
    outR(i) = makeWeightedCorr1(Q{h(i)}, reshape(Pr(bID == outID(i), :), [], 1));
end

    
end

function out = makeWeightedCorr1(xy, w);
%function out = makeWeightedCorr1(xy, w);

mxy = sum(xy.*repmat(w, 1, 2)/sum(w));
covxy = sum(w.*(xy(:, 1) - mxy(1)).*(xy(:, 2) - mxy(2)))/sum(w);
covxx = sum(w.*(xy(:, 1) - mxy(1)).*(xy(:, 1) - mxy(1)))/sum(w);
covyy = sum(w.*(xy(:, 2) - mxy(2)).*(xy(:, 2) - mxy(2)))/sum(w);
out = covxy/(sqrt(covyy*covxx));
end
