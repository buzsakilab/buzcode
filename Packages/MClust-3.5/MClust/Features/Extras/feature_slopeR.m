function [slopeData, slopeNames, slopePars] = feature_SlopeR(V, ttChannelValidity, Params)

% MClust% [slopeData, slopeNames, slopePars] = feature_SlopeR(V, ttChannelValidity, Params)% Calculate maximum rising slope feature. Performs cubic spline
% interpolation to improve resolution.%% INPUTS%    V = TT tsd%    ttChannelValidity = nCh x 1 of booleans%% OUTPUTS%    slopeData - nSpikes x nCh valley values%    slopeNames - "valley1to10: Ch"
%
% JCJ Sept 2002%

TTData = Data(V);

[nSpikes, nCh, nSamp] = size(TTData);

f = find(ttChannelValidity);

slopeData = zeros(nSpikes, length(f));

slopeNames = cell(length(f), 1);
slopePars = {};
slopeData = squeeze(max(diff(TTData(:, f, :),1,3), [], 3));

for iCh = 1:length(f)
   slopeNames{iCh} = ['SlopeR: ' num2str(f(iCh))];
end
