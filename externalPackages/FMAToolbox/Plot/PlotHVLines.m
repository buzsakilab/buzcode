function p = PlotHVLines(positions,direction,varargin)

%PlotHVLines - Plot vertical (resp. horizontal) lines at listed x (resp. y).
%
%  USAGE
%
%    p = PlotHVLines(positions,direction,options)
%
%    positions      list of abscissae/ordinates
%    direction      optional direction: 'h' or 'v' (default = 'v')
%    <options>      options for function <a href="matlab:help plot">plot</a>
%

% Copyright (C) 2008-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
 	error('Incorrect number of parameters (type ''help <a href="matlab:help PlotHVLines">PlotHVLines</a>'' for details).');
end
if nargin < 2,
	direction = 'v';
else
	direction = lower(direction);
end
if min(size(positions)) > 2,
 	error('List of abscissae/ordinates is not a vector (type ''help <a href="matlab:help PlotHVLines">PlotHVLines</a>'' for details).');
else
	positions = positions(:);
end

hold on;
if strcmp(direction,'v'),
	yLim = ylim;
	for i = 1:size(positions,1),
		plot([positions(i,1) positions(i,1)],yLim,varargin{:});
	end
else
	xLim = xlim;
	for i = 1:size(positions,1),
		plot(xLim,[positions(i,1) positions(i,1)],varargin{:});
	end
end