function [X,Y,V,GoodRanges] = LoadPosition_Wrapper(fbasename, scaleFactor)


% loads position and speed from a position file file (ending in .whl)
%
% USAGE
%     [X,Y,V,GoodRanges] = LoadPosition_Wrapper(fbasename)
%     
% INPUT:
%     fbasename: session file basename
%     scaleFactor: scaling factor to convert pixels to cm (in units of cm/pixel)
%	
% OUTPUT:
%     X: a tsd object of x position values
%     Y: a tsd object of y position values
%     V: a tsd object of velocity values
%     GoodRanges: a intervalSet object where LEDs were successfully detected

% Adrien Peyrache, 2011
% edited by Luke Sjulson to add scaleFactor and use units of seconds, 2015-11


smoothWidth = 30;

[whl,t,GoodRanges] = LoadPosition(fbasename);
if size(whl)>2
    X = nanmean([whl(:,1),whl(:,3)]')';
    Y = nanmean([whl(:,2),whl(:,4)]')';
else
    X = whl(:,1);
    Y = whl(:,2);
end

% convert to units of cm
if nargin<2, scaleFactor = 1; end
X = X.*scaleFactor;
Y = Y.*scaleFactor;

dx = diff(X);
dy = diff(Y);
dt = median(diff(t));
v = [sqrt(dx.^2+dy.^2);0];

gw = gausswin(smoothWidth);
gw = gw/sum(gw(:));
warning off
goodEp = intervalSet(GoodRanges(:,1),GoodRanges(:,2));
goodEp = mergeCloseIntervals(goodEp,1);
warning on

v(isnan(v))=0;
V = convn(v,gw,'same');
V = V/dt;

X = tsd(t,X);
Y = tsd(t,Y);
V = tsd(t, V);

GoodRanges = intervalSet(GoodRanges(:,1),GoodRanges(:,2));
