function [h X Y] = hist2(Data, XBins, YBins)
% out = hist2(Data, XBins, YBins)
%
% Makess a 2d histogram of the data.
% XBins and YBins are optional arguments
% which give the number of grid segments.
% - default is 50.

if length(XBins) == 1
    MinX = min(Data(:,1));
    MaxX = max(Data(:,1));
else
    MinX = XBins(1);
    MaxX = XBins(end);
    XBins = length(XBins);
end
if length(YBins) == 1
   MinY = min(Data(:,2));
   MaxY = max(Data(:,2));
else
    MinY = YBins(1);
    MaxY = YBins(end);
    YBins = length(YBins);
end

XBin = floor(1 + XBins*(Data(:,1) - MinX) / (MaxX - MinX));
YBin = floor(1 + YBins*(Data(:,2) - MinY) / (MaxY - MinY));

ix = isnan(XBin) | isnan(YBin);
XBin(ix) = [];
YBin(ix) = [];

h = zeros(XBins, YBins);

XBin(find(XBin == XBins+1)) = XBins;
YBin(find(YBin == YBins+1)) = YBins;

for i = 1:size(XBin,1)
    h(XBin(i), YBin(i)) = h(XBin(i), YBin(i)) + 1;
end;

X = [MinX:(MaxX-MinX)/(XBins-1):1.000000001*MaxX];
Y = [MinY:(MaxY-MinY)/(YBins-1):1.000000001*MaxY];


%  imagesc(h)

% Written by Kenneth D. Harris, adapted by Adrien Peyrache
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu