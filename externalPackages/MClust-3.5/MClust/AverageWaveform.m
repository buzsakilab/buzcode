function [mWV, sWV] = AverageWaveform(WV,varargin)

% AverageWaveform(wv)
%
%
% INPUTS
%    wv -- tsd of tt or waveform data
%
% OUTPUTS
%    mWV - mean waveform 4 x 32
%    sWV - stddev waveform 4 x 32
%
%
% ADR 1998
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.
%
% modified by ncst 23 Apr 03 to be compatible with n samples 
% (instead of 32)

DoPlotYN = 'no';

extract_varargin;

WVD = Data(WV);

nWVSamples = size(WVD,3);

mWV = zeros(4,nWVSamples);
sWV = zeros(4,nWVSamples);

mWV = squeeze(mean(WVD,1));
sWV = squeeze(std(WVD,1));

if nargout == 0 || strcmp(DoPlotYN,'yes')  % no output args, plot it
    if strcmp(DoPlotYN,'yes')
        AveWVFig = figure;
    end
    for it = 1:4
        xrange = ((nWVSamples + 2) * (it-1)) + (1:nWVSamples); 
        if strcmp(DoPlotYN,'yes')
            figure(AveWVFig);
        end
        hold on;
        plot(xrange, mWV(it,:));
        errorbar(xrange,mWV(it,:),sWV(it,:)); 
    end
    axis off
    axis([0 4*(nWVSamples + 2)])
    title('Average Waveform');
    hold off
end