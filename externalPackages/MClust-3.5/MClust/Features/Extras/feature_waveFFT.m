function [FFTData, FFTNames, FFTPars] = feature_FFT(V, ttChannelValidity, Param)

% MClust
% [Data, Names] = feature_FFT(V, ttChannelValidity)
% Calculate FFT feature value for each channel
% 
%
% INPUTS
%    V = TT tsd
%    ttChannelValidity = nCh x 1 of booleans
%
% OUTPUTS
%    Data - nSpikes x nCh of weighted sum of fast fourier transform of each spike
%    Names - "waveFFT: Ch"

% ncst Jan 2002%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

TTData = Data(V);
[nSpikes, nCh, nSamp] = size(TTData);
f = find(ttChannelValidity);
FFTData = zeros(nSpikes, length(f));
FFTNames = cell(length(f), 1);
FFTPars = {};

for iCh = 1:length(f)
    Y = fft(squeeze(TTData(:,f(iCh),:))',32);
    Pyy = Y.*conj(Y)/32;
    WeightMatrix = repmat(([1:16 16:-1:1])',1,length(Pyy(1,:)));
    SumPyy = sum(Pyy);
    FFTData(:, iCh) = (sum(Pyy.*WeightMatrix)./SumPyy)';     FFTNames{iCh} = ['waveFFT: ' num2str(f(iCh))];
end
