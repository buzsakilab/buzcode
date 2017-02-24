function PlotColorCurves(curves,range,varargin)

%PlotColorCurves - Plot a series of curves as an aggregated color map.
%
%  Plot a series of curves as an aggregated color map. In addition, one or more
%  reference x locations can be represented as vertical lines: examples include
%  time zero for cross-correlograms, baseline frequency for power spectra, etc.
%  Finally, colored markers can mark locations of interest on each individual
%  curve: examples include local minima, zero crossings, etc.
%
%  USAGE
%
%    PlotColorCurves(curves,range,x1,style1,x2,style2,...,<options>)
%
%    curves         MxN matrix consisting of M curves (N bins each)
%    range          optional x values for first and last bin centers
%                   (if subsequent parameters are provided, this is
%                   no longer optional)
%    x1,x2...       optional locations of interest (see below)
%    style1...      optional drawing style, e.g. 'w--', 'k.', etc.
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'cutoffs'     lower and upper cutoff values ([] = autoscale, default)
%     'hgamma'      gamma-like correction for hue (1 = no correction, default)
%     'colorspace'  'HSV' or 'RGB' (default = 'RGB')
%     'bar'         draw a color bar (default = 'off')
%     'type'        either 'linear' or 'circular' (default 'linear')
%    =========================================================================
%
%  NOTE
%
%    Locations of interest can be either scalar values (for vertical lines) or
%    vectors of length M (for individual markers).
%

% Copyright (C) 2012-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

colors = 'kbgrw';

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotColorCurves">PlotColorCurves</a>'' for details).');
end

% Check parameters
if ~isdmatrix(curves),
	error('Incorrect curves (type ''help <a href="matlab:help PlotColorCurves">PlotColorCurves</a>'' for details).');
end
if nargin < 2,
	range = [0 size(curves,2)];
elseif ~isdvector(range,'#2','<'),
	error('Incorrect range (type ''help <a href="matlab:help PlotColorCurves">PlotColorCurves</a>'' for details).');
end

options = {};
last = length(varargin);
% Parse parameter list: find options for PlotColorMap
for i = 1:length(varargin),
	if ischar(varargin{i}),
		switch(lower(varargin{i})),
			case {'cutoffs','gamma','hgamma','colorspace','bar','type'},
				options = {varargin{i:end}};
				last = i-1;
				break;
		end
	end
end

% Plot color curves
[m,n] = size(curves);
if isempty(options),
	PlotColorMap(curves,'x',linspace(range(1),range(2),n));
else
	PlotColorMap(curves,'x',linspace(range(1),range(2),n),options{:});
end
hold on;

i = 1;
c = 1;
while i <= last,
	x = varargin{i};
	if isdscalar(x),
		if i+1 <= length(varargin) && ischar(varargin{i+1}),
			style = varargin{i+1};
			PlotHVLines(x,'v',style);
			i = i + 2;
		else
			PlotHVLines(x,'v','w--');		
			i = i + 1;
		end
	elseif isdvector(x,['#' int2str(m)]),
		if i+1 <= length(varargin) && ischar(varargin{i+1}),
			style = varargin{i+1};
			plot(x,1:m,style,'markersize',8);
			i = i + 2;
		else
			plot(x,1:m,['.' colors(k)],'markersize',8);
			i = i + 1;			
			k = k + 1;
			if k > length(colors), k = 1; end
		end
	else
		error('Incorrect locations (type ''help <a href="matlab:help PlotColorCurves">PlotColorCurves</a>'' for details).');	
	end
end
