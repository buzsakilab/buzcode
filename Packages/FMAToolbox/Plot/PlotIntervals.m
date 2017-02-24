function PlotIntervals(intervals,style,direction)

%PlotIntervals - Plot vertical bars or rectangles to show interval limits.
%
% Given a list of intervals [start stop], draw a green vertical bar at
% the beginning of each interval, and a red vertical bar at the end of
% each interval or a grey rectangle representing the interval.
%
%  USAGE
%
%    PlotIntervals(intervals,style,direction)
%
%    intervals      list of [start stop] intervals
%    style          optional style: 'bars' or 'rectangles'
%    direction      optional direction: 'h' or 'v' (default = 'v')

% Copyright (C) 2008-2011 by Gabrielle Girardeau & Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
end
if size(intervals,2) ~= 2,
  error('Incorrect list of intervals (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
end
if nargin < 2,
	style = 'rectangles';
else
	style = lower(style);
end
if nargin < 3,
	direction = 'v';
else
	direction = lower(direction);
end

if ~isstring_FMAT(style,'bars','rectangles'),
	error('Incorrect value for property ''style'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
end
if ~isstring_FMAT(direction,'h','v'),
	error('Incorrect value for property ''direction'' (type ''help <a href="matlab:help PlotIntervals">PlotIntervals</a>'' for details).');
end

hold on;
yLim = ylim;
xLim = xlim;
if strcmp(style,'bars'),
	for i = 1:size(intervals,1),
		if strcmp(direction,'v'),
			plot([intervals(i,1) intervals(i,1)],yLim,'Color',[0 0.75 0]);
			plot([intervals(i,2) intervals(i,2)],yLim,'Color',[0.9 0 0]);
		else
			plot(xLim,[intervals(i,1) intervals(i,1)],'Color',[0 0.75 0]);
			plot(xLim,[intervals(i,2) intervals(i,2)],'Color',[0.9 0 0]);
		end
	end
else
	for i=1:size(intervals,1),
		if strcmp(direction,'v'),
			rec = rectangle('Position',[intervals(i,1) yLim(1) intervals(i,2)-intervals(i,1) yLim(2)-yLim(1)],'FaceColor',[0.9 0.9 0.9],'LineStyle','none');
		else
			rec = rectangle('Position',[xLim(1) intervals(i,1) xLim(2)-xLim(1) intervals(i,2)-intervals(i,1)],'FaceColor',[0.9 0.9 0.9],'LineStyle','none');
		end
		uistack(rec,'bottom');
	end
end


