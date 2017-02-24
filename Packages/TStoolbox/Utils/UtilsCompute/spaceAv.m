function [h X Y] = hist2(Pos, tPos, tsa, XBins, YBins)

% USAGE :
% [out binX binY] = hist2(Data, XBins, YBins)
%
% Makess a 2d histogram of the data.
%  INPUT:
%  	Pos: a n-by-2 matrix of postions
%  	tsa: a tsd to be space averaged
%  	XBins and YBins are optional arguments
%	   which give the number of grid segments.
% 	   - default is 50.
%  OUTPUT:
%  	h: the space average matrix
%  	binX, binY : the vectors which each elements are the centers of the space bins
%  
% Written by Adrien Peyrache from hist2 written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html


if (nargin<2) XBins = 50; end;
if (nargin<3) YBins = 50; end;

MinX = min(Pos(:,1));
MaxX = max(Pos(:,1));
MinY = min(Pos(:,2));
MaxY = max(Pos(:,2));

binPos = zeros(XBins, YBins);
h = zeros(XBins, YBins);

d = Data(Restrict(tsa,ts(tPos)));

XBin = floor(1 + XBins*(Pos(:,1) - MinX) / (MaxX - MinX));
YBin = floor(1 + YBins*(Pos(:,2) - MinY) / (MaxY - MinY));

XBin(find(XBin == XBins+1)) = XBins;
YBin(find(YBin == YBins+1)) = YBins;

for i = 1:size(Pos,1)
	binPos(XBin(i), YBin(i)) = binPos(XBin(i), YBin(i))+1;
end;

for i = 1:size(Pos,1)
	h(XBin(i), YBin(i)) = h(XBin(i), YBin(i)) + d(i);
end;

warning off

	h= h./binPos;

warning on

X = [MinX:(MaxX-MinX)/(XBins-1):MaxX];
Y = [MinY:(MaxY-MinY)/(YBins-1):MaxY];

%  imagesc(h)
