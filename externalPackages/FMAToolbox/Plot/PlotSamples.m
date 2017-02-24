function p = PlotSamples(X,varargin)

%PlotSamples - Plot samples (multiple time series).
%
% Plot columns 2...M of a matrix as a function of its first column.
%
%  USAGE
%
%    p = PlotSamples(X,<options>,<options2>)
%
%    X              list of (t,x1,x2...xN) data to plot (N variables)
%    <options>      optional list of property-value pairs (see table below)
%    <options2>     options for function <a href="matlab:help plot">plot</a>
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'spacing'     vertical spacing between variables (default = 1 for point
%                   processes, 0 for continuous data)
%    =========================================================================
%

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotSamples">PlotSamples</a>'' for details).');
end

% Default values
pointProcess = false;
spacing = [];
v = {};

% Parse parameter list
n = 1;
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotSamples">PlotSamples</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'spacing',
			spacing = varargin{i+1};
			if ~isdscalar(spacing),
				error('Incorrect value for property ''spacing'' (type ''help <a href="matlab:help PlotSamples">PlotSamples</a>'' for details).');
			end

		case 'type',
			type = lower(varargin{i+1});
			if ~isstring_FMAT(type,'continuous','point'),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help PlotSamples">PlotSamples</a>'' for details).');
			end
			pointProcess = strcmp(type,'point');

		otherwise,
			v = varargin{i:end};
			if ~isa(v,'cell'), v = {v}; end
			break;
	end
end

if isempty(spacing),
	if pointProcess,
		spacing = 1;
	else
		spacing = 0;
	end
end

% Plot
hold on;
if pointProcess,
	for i = 1:size(X,2),
		PlotTicks(X(:,i),(i-1)*spacing,'direction','v',v{:});
	end
else
	for i = 2:size(X,2),
		plot(X(:,1),(i-1)*spacing+X(:,i),v{:});
	end
end

