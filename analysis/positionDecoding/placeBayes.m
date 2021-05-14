function [Pr, prMax] = placeBayes(Cr, rateMap, binLength)
%function out = placeBayes1(Cr, rateMap, binLength)
% Inputs:
% Cr = [nTemporalBin X nCell] matrix of binned firing rates
% rateMap = [nSpatialBin X nCell] firing rate 'template'
% binLength = scalar of the duration of the bins in 'Cr'
% Outputs:
% Pr = [nTemporalBin X nSpatialBins] matrix of posterior probabilities
% prMax = the spatial bin with higher spatial probabilities for each
%   temporalBin in Cr (note, ties go to the lower numbered bins as 
%   consistent with the behavior of the second output of built in function
%   'max')


Cr = Cr*binLength;
rateMap = rateMap';
term2 = exp((-1)*binLength*sum(rateMap'));
mp = 1/size(rateMap, 2);
Pr = [];

c =  repmat(Cr, [1, 1, size(rateMap, 1)]);
b = repmat(rateMap', [1, 1, size(c, 1)]);
b = shiftdim(b, 2);


u = mp*prod(b.^c, 2);
u = squeeze(u);
Pr = u.*repmat(term2, size(u, 1), 1);
Pr = Pr./repmat(sum(Pr, 2), 1, size(Pr, 2));

[~, m] = max(Pr');
prMax = m';

if sum(sum(isinf(Pr))) > 0
    error('Do Not Approach the Infitnite');
end

% if sum(sum(isnan(Pr))) > 0
%     error('What is ''not a nubmer''?');
% end
