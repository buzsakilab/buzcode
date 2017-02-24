function [energyData, energyNames, energyPars] = feature_EnergyD1(V, ttChannelValidity, Params)

% MClust
% [Data, Names, Params] = feature_EnergyD1(V, ttChannelValidity, Params)
% Calculate energy feature for the first derivative of each channel.
% Normalizes for number of samples in waveform.
%
% INPUTS
%    V = TT tsd
%    ttChannelValidity = nCh x 1 of booleans%    Params   = feature paramters  (none for energy) %% OUTPUTS%    energyData - nSpikes x nCh of energy INSIDE curve (below peak and above valley) of each spike%    energyNames - "energy: Ch"

% JCJ Nov 2003

TTData = Data(V);

[nSpikes, nCh, nSamp] = size(TTData);

f = find(ttChannelValidity);

energyData = zeros(nSpikes, length(f));

energyNames = cell(length(f), 1);
energyPars = {};
energyData = squeeze(sqrt(sum(diff(TTData(:, f, :),1, 3).^2,3)))./sqrt(nSamp);

for iCh = 1:length(f)
   energyNames{iCh} = ['EnergyD1: ' num2str(f(iCh))];
end
