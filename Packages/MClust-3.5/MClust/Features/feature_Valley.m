function [valleyData, valleyNames, valleyPars] = feature_Valley(V, ttChannelValidity, Params)

% MClust% [valleyData, valleyNames] = feature_Valley(V, ttChannelValidity)% Calculate valley feature max value for each channel%% INPUTS%    V = TT tsd%    ttChannelValidity = nCh x 1 of booleans%% OUTPUTS%    Data - nSpikes x nCh valley values%    Names - "valley: Ch"
%
% ADR April 1998%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.


TTData = Data(V);
[nSpikes, nCh, nSamp] = size(TTData);

f = find(ttChannelValidity);

valleyData = zeros(nSpikes, length(f));
valleyNames = cell(length(f), 1);
valleyPars = {};

for iCh = 1:length(f)
   valleyData(:, iCh) = min(squeeze(TTData(:, f(iCh), :)), [], 2);
   valleyNames{iCh} = ['Valley: ' num2str(f(iCh))];
end
