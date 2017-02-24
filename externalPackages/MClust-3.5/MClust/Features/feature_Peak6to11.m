function [PeakData, PeakNames,PeakPars] = feature_PEAK6to11(V, ttChannelValidity, Params)

% MClust
% [PeakData, PeakNames] = feature_PEAK6to11(V, ttChannelValidity)
% Calculate peak feature max value for each channel
%
% INPUTS
%    V = TT tsd
%    ttChannelValidity = nCh x 1 of booleans
%
% OUTPUTS
%    Data - nSpikes x nCh peak values
%    Names - "Peak6to11: Ch"
%
% ADR April 1998
% version M1.0
% RELEASED as part of MClust 2.0
% See standard disclaimer in Contents.m

TTData = Data(V);
[nSpikes, nCh, nSamp] = size(TTData);

f = find(ttChannelValidity);

PeakData = zeros(nSpikes, length(f));
PeakNames = cell(length(f), 1);
PeakPars = {};
PeakData = squeeze(max(TTData(:, f, 6:11), [], 3));

for iCh = 1:length(f)
   PeakNames{iCh} = ['Peak6to11: ' num2str(f(iCh))];
end
