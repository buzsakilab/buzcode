function [PeakData, PeakNames, PeakPars] = feature_PeakValleyDiff(V, ttChannelValidity,Params)

% MClust
% [PvData, PvNames] = feature_PeakValleyDiff(V, ttChannelValidity)
% Calculate pv feature max value for each channel
%
% INPUTS
%    V = TT tsd
%    ttChannelValidity = nCh x 1 of booleans
%
% OUTPUTS
%    Data - nSpikes x nCh pv values
%    Names - "Pv: Ch"
%
% ADR April 1998
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

TTData = Data(V);

[nSpikes, nCh, nSamp] = size(TTData);

f = find(ttChannelValidity);

PeakData = zeros(nSpikes, length(f));

PeakNames = cell(length(f), 1);
PeakPars = {};
PeakData = squeeze(max(TTData(:, f, :), [], 3) - min(TTData(:, f, :), [], 3));

for iCh = 1:length(f)
   PeakNames{iCh} = ['PeakValleyDiff: ' num2str(f(iCh))];
end
