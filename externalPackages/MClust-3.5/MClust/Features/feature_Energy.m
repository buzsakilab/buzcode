function [energyData, energyNames, energyPars] = feature_Energy(V, ttChannelValidity, Params)

% MClust
% [Data, Names, Params] = feature_Energy(V, ttChannelValidity, Params)
% Calculate energy feature max value for each channel. Normalizes for #
% samples in waveform.
%
% INPUTS
%    V = TT tsd
%    ttChannelValidity = nCh x 1 of booleans%    Params   = feature paramters  (none for energy) %% OUTPUTS%    Data - nSpikes x nCh of energy INSIDE curve (below peak and above valley) of each spike%    Names - "energy: Ch"

% ADR April 1998
% JCJ Nov 2003 Modified to use cubic spline


TTData = Data(V);

[nSpikes, nCh, nSamp] = size(TTData);

f = find(ttChannelValidity);

energyData = zeros(nSpikes, length(f));

energyNames = cell(length(f), 1);
energyPars = {};
energyData = squeeze(sqrt(sum(TTData(:, f, :).^2,3)))./sqrt(nSamp);

for iCh = 1:length(f)
   energyNames{iCh} = ['Energy: ' num2str(f(iCh))];
end
