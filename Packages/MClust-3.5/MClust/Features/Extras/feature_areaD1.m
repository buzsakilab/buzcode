function [AreaData, AreaNames, AreaPars] = feature_AreaD1(V, ttChannelValidity, Params)

% MClust% [Data, Names] = feature_AreaD1(V, ttChannelValidity)% Calculate area feature for 1st derivative for each channel. Normalizes
% for number of samples in waveform.%% INPUTS%    V = TT tsd%    ttChannelValidity = nCh x 1 of booleans%% OUTPUTS%    Data - nSpikes x nCh of area INSIDE curve (below peak and above valley) of each spike%    Names - "Area: Ch"
%
% JCJ Nov 2003%

TTData = Data(V);

[nSpikes, nCh, nSamp] = size(TTData);

f = find(ttChannelValidity);

AreaData = zeros(nSpikes, length(f));

AreaNames = cell(length(f), 1);
AreaPars = {};
AreaData = squeeze(sum(diff(TTData(:, f, :),1, 3),3))./nSamp;

for iCh = 1:length(f)
   AreaNames{iCh} = ['AreaD1: ' num2str(f(iCh))];
end
