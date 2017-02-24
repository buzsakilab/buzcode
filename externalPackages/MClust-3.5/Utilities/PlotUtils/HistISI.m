function [H, binsUsed] = HistISI(TS, varargin)

% H = HistISI(TS, parameters)
%  H = HistISI(TS, 'maxLogISI','maxLogISI',5)      for fixed upper limit 10^5 msec (or 100 sec)
%
% INPUTS:
%      TS = a single ts object
%
% OUTPUTS:
%      H = histogram of ISI
%      N = bin centers
%
% PARAMETERS:
%     nBins 500
%     maxLogISI variable
%     minLogISI 
%
% If no outputs are given, then plots the figure directly.
%
% Assumes TS is in seconds or timestamps!
%
% ADR 1998
% version L5.3
% RELEASED as part of MClust 2.0
% See standard disclaimer in Contents.m
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.
%
% ADR 2 DEC 2003 fixed capitalizations

%--------------------
if ~isa(TS, 'ts'); 
    error('Type Error: input is not a ts object.'); 
end
epsilon = 1e-100;
nBins = 500;
maxLogISI = 0;
minLogISI = [];
DoPlotYN = 'no';

extract_varargin;

%--------------------
% Assumes data is passed in in either seconds or timestamps
% if the maximum of TS is less than 10^6, then it assumes we need
% to convert to timestamps
%
% modified ncst 13 Apr 03, to use 1e-6 rather than 1e-7

if max(Data(TS)) < 1e6
    TS = ts(Data(TS)*10000);
end

ISI = diff(Data(TS)/10) + epsilon;
if length(ISI) == 0
   warning('ISI contains no data!');
   ISI = 1;
end   
if ~isreal(log10(ISI))
   warning('ISI contains negative differences; log10(ISI) is complex.');
   complexISIs = 1;
else
   complexISIs = 0;
end

if maxLogISI == 0
    maxLogISI = max(real(log10(ISI)))+1;     
end

if isempty(minLogISI)
    minLogISI = floor(min(real(log10(ISI))));
end

H = ndhist(log10(ISI)', nBins, minLogISI, maxLogISI);
binsUsed = logspace(minLogISI,maxLogISI,nBins);

%-------------------
if nargout == 0 | strcmp(DoPlotYN,'yes')
    if strcmp(DoPlotYN,'yes')
        ISIFig = figure;
    end
    
    if strcmp(DoPlotYN,'yes')
        figure(ISIFig);
    end
    plot(binsUsed, H);
    if complexISIs
        xlabel('ISI, ms.  WARNING: contains negative components.');
    else
        xlabel('ISI, ms');
    end
    set(gca, 'XScale', 'log', 'XLim', [10^minLogISI 10^maxLogISI]);
    set(gca, 'YTick', max(H));
    
    % draw line at ms
    if strcmp(DoPlotYN,'yes')
        figure(ISIFig);
    end
    hold on
    plot([1 1], get(gca, 'YLim'), 'r:')
    hold off
end

   


   