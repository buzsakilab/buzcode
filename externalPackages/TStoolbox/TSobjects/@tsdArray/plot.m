function h = plot(S,tickson)
% Plots a tsdArray as raster (based on time only, no values used).
% Brendon Watson 2015

if ~exist('tickson','var')
    tickson = 0;
end

x = [];
y = [];
S = cellArray(S);

for a = 1:length(S);
    tx = Range(S{a},'s');
    ty  = a*ones(size(tx));
    x = cat(1,x,tx);
    y = cat(1,y,ty);
end

if tickson
    y1 = y-.4;
    y2 = y+.4;
    plot([x x]',[y1 y2]')
else
    plot(x,y,'.')
    axis tight
end

yl = ylim;
ylim([yl(1)-.5 yl(2)+.5])
