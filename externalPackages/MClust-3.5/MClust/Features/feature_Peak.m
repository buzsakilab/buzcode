function [PeakData, PeakNames,PeakPars] = feature_Peak(V, ttChannelValidity, Params)

% MClust
% [PeakData, PeakNames] = feature_Peak(V, ttChannelValidity)
% Calculate peak feature max value for each channel
%
% INPUTS
%    V = TT tsd
%    ttChannelValidity = nCh x 1 of booleans
%
% OUTPUTS
%    Data - nSpikes x nCh peak values
%    Names - "Peak: Ch"
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
for iCh = 1:length(f)
	PeakData(:,iCh) = squeeze(max(TTData(:, f(iCh), :), [], 3));
	PeakNames{iCh} = ['_Peak: ' num2str(f(iCh))];
end
