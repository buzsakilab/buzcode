function plot(is,col)
% makes a plot to represent an intervalset but simply plotting horizontal
% lines based on the start and stop time of each interval.  Units of plot
% is in seconds.
%
% INPUTS:
% is - an intervalSet
% col - optional color input for input to native "plot" command
%
% OUTPUTS:
% (only graphical)

% Brendon Watson 2014

% SampFreq = 10000;

spans = [Start(is,'s') End(is,'s')];
% spans = spans/SampFreq;

if exist('col','var')
%     plot(spans',zeros(size(spans')),'color',col)
    plot(spans',zeros(size(spans')),'color',col,'linewidth',3)
else
%     plot(spans',zeros(size(spans')))
    plot(spans',zeros(size(spans')),'linewidth',3)
end