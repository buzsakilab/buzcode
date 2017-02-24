function h = PlotRepeat(X,Y,pattern,varargin)

%PlotRepeat - Plot repeated (tiled) copies of the data.
%
%  USAGE
%
%    h = PlotRepeat(X,Y,pattern,<options>,<options2>)
%
%    X              abscissae
%    Y              ordinates
%    pattern        'xxy','xyy' or 'xxyy' depending on repeated axes
%    <options>      optional list of property-value pairs (see table below)
%    <options2>     additional options for 'plot' (for 'curve' mode only)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'mode'        plot type, 'curve' (default) or 'bar' (only for xxy)
%     'xlim'        [min max] for X
%     'ylim'        [min max] for Y
%    =========================================================================

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
xLims = [];
yLims = [];
parameters = {};
mode = 'curve';

% Check number of parameters
if nargin < 3,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotRepeat">PlotRepeat</a>'' for details).');
end

% Check pattern
if ~isstring_FMAT(pattern,'xxy','xyy','xxyy'),
  error('Incorrect pattern (type ''help <a href="matlab:help PlotRepeat">PlotRepeat</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotRepeat">PlotRepeat</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'mode',
			mode = varargin{i+1};
			if ~isstring_FMAT(mode,'curve','bar'),
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help PlotRepeat">PlotRepeat</a>'' for details).');
			end
		case 'xlim',
			xLims = varargin{i+1};
			if ~isdvector(xLims,'<','#2'),
				error('Incorrect value for property ''xlim'' (type ''help <a href="matlab:help PlotRepeat">PlotRepeat</a>'' for details).');
			end
		case 'ylim',
			yLims = varargin{i+1};
			if ~isdvector(yLims,'<','#2'),
				error('Incorrect value for property ''ylim'' (type ''help <a href="matlab:help PlotRepeat">PlotRepeat</a>'' for details).');
			end
		otherwise,
			parameters = varargin{i:end};
			if ~isa(parameters,'cell'), parameters = {parameters}; end
			break;
	end
end

% Check pattern and mode consistency
if strcmp(mode,'bar') && ~strcmp(pattern,'xxy'),
  error('Bar plots can only be used for xxy pattern (type ''help <a href="matlab:help PlotRepeat">PlotRepeat</a>'' for details).');
end

X = X(:);
Y = Y(:);

if isempty(X) || isempty(Y),
	return;
end

if isempty(xLims),
	xLims = [min(X) max(X)];
end
if isempty(yLims),
	yLims = [min(Y) max(Y)];
end

hold on;
switch pattern,
	case 'xxy',
		if strcmp(mode,'curve'),
			h = plot([X;X+diff(xLims)],[Y;Y],parameters{:});
		else
			[X,i] = unique([X;X+diff(xLims)]);
			Y = [Y;Y];
			Y = Y(i);
			h = bar(X,Y,'hist');
		end
		xlim(ConsolidateIntervals([xlim;xLims;xLims+diff(xLims)]));
		ylim(ConsolidateIntervals([ylim;yLims]));
	case 'xyy',
		h = plot([X;X],[Y;Y+diff(yLims)],parameters{:});
		xlim(ConsolidateIntervals([xlim;xLims]));
		ylim(ConsolidateIntervals([ylim;yLims;yLims+diff(yLims)]));
	case 'xxyy',
		h = plot([X;X+diff(xLims);X;X+diff(xLims)],[Y;Y;Y+diff(yLims);Y+diff(yLims)],parameters{:});
		xlim(ConsolidateIntervals([xlim;xLims;xLims+diff(xLims)]));
		ylim(ConsolidateIntervals([ylim;yLims;yLims+diff(yLims)]));
end
