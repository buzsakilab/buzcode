function [map,stats] = Map(v,z,varargin)

%Map - Map z on (x,y) where x, y and z are time-varying variables (samples).
%
%  Compute a continuous map, where one time-varying variable z is represented
%  as a function of one or two time-varying variables x and y. The variable z
%  can either be a point process (typically, a list of spike timestamps) or a
%  continuous measure (e.g. the instantaneous velocity of the animal, the
%  spectral power of an LFP channel in a given frequency band, the coherence
%  between two oscillating LFP channels, etc.) Typical examples of x and y
%  include spatial coordinates and angular directions.
%
%  An occupancy map is also computed.
%
%  USAGE
%
%    map = Map([t x y],z,<options>)
%
%    t              timestamps for x and y
%    x              x values in [0,1]
%    y              optional y values in [0,1]
%    z              list of timestamps or (timestamp,value) pairs
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'       number of horizontal and vertical bins (default = [50 50])
%     'minTime'     minimum time spent in each bin (in s, default = 0)
%     'maxGap'      z values recorded during time gaps between successive (x,y)
%                   samples exceeding this threshold (e.g. undetects) will not
%                   be interpolated; also, such long gaps in (x,y) sampling
%                   will be clipped to 'maxGap' to compute the occupancy map
%                   (default = 0.100 s)
%     'type'        'linear' if z is linear (default), 'circular' otherwise
%    =========================================================================
%
%  OUTPUT
%
%    map.x          x bins
%    map.y          y bins
%    map.z          average map (z continuous)
%                   or rate map (z point process)
%    map.count      count map (z point process)
%    map.time       occupancy map (in s)
%
%  NOTES
%
%    x values are arranged in columns and y values in rows in all output matrices
%    (e.g. 'map.z').
%
%  SEE
%
%    See also MapStats, FiringMap, PlotColorMap, Accumulate.

% Copyright (C) 2002-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Map">Map</a>'' for details).');
end

% Check parameter sizes
if size(v,2) < 2,
	error('Parameter ''[t x y]'' should have at least 2 columns (type ''help <a href="matlab:help Map">Map</a>'' for details).');
end
if (size(z,2) < 1 || size(z,2) > 2) && ~isempty(z),
	error('Parameter ''z'' should have 1 or 2 columns (type ''help <a href="matlab:help Map">Map</a>'' for details).');
end

% Default values
maxGap = 0.1;
map.x = [];
map.y = [];
map.count = [];
map.time = [];
map.z = [];
smooth = 2;
nBins = 50;
minTime = 0;
type = 'c';

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Map">Map</a>'' for details).']);
	end
	switch(lower(varargin{i})),

		case 'smooth',
			smooth = varargin{i+1};
			if ~isdvector(smooth,'>=0') | length(smooth) > 2,
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
			end

		case 'nbins',
			nBins = varargin{i+1};
			if ~isivector(nBins,'>0') | length(nBins) > 2,
				error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
			end

		case 'mintime',
			minTime = varargin{i+1};
			if ~isdscalar(minTime,'>=0'),
				error('Incorrect value for property ''minTime'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
			end

		case 'maxgap',
			maxGap = varargin{i+1};
			if ~isdscalar(maxGap,'>=0'),
			error('Incorrect value for property ''maxGap'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
			end

		case 'type',
			type = lower(varargin{i+1});
			if ~isstring_FMAT(type,'circular','linear'),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help Map">Map</a>'' for details).');
			end

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Map">Map</a>'' for details).']);

  end
end

if isempty(v), return; end

% Some info about x, y and z
pointProcess = (isempty(z) | size(z,2) == 1);
t = v(:,1);
x = v(:,2);
if size(v,2) >= 3,
	y = v(:,3);
else
	y = [];
end

% Make sure x and y are normalized
if max(x) > 1 || min(x) < 0,
	x = ZeroToOne(x);
	warning('Parameter ''x'' should contain values in [0 1]. The data will now be transformed accordingly.');
end
if ~isempty(y),
	if max(y) > 1 || min(y) < 0,
		y = ZeroToOne(y);
		warning('Parameter ''y'' should contain values in [0 1]. The data will now be transformed accordingly.');
	end
end

% Number of bins for x and y
nBinsX = nBins(1);
if length(nBins) == 1,
	nBinsY = nBinsX;
	nBins(2) = nBins;
else
	nBinsY = nBins(2);
end

% Bin x and y
x = Bin(x,[0 1],nBinsX);
if ~isempty(y),
	y = Bin(y,[0 1],nBinsY);
end

% Duration for each (X,Y) sample (clipped to maxGap)
dt = diff(t);dt(end+1)=dt(end);dt(dt>maxGap) = maxGap;

if pointProcess,
	% Count occurrences for each (x,y) timestamp
	n = CountInIntervals(z,[t t+dt]);
%  	n = CountInIntervals(z,[t(1:end-1) t(2:end)]);
%  	n(end+1) = 0;
else
	% Interpolate z at (x,y) timestamps
	[z,discarded] = Interpolate(z,t,'maxGap',maxGap);
	if isempty(z), return; end
	if strcmp(type,'circular'),
		range = isradians(z(:,2));
		z(:,2) = exp(j*z(:,2));
	end
	n = 1;
end

% Computations
if isempty(y),
	% 1D (only x)
	map.x = linspace(0,1,nBinsX);
	map.count = Smooth(Accumulate(x,n,nBinsX),smooth)';
	map.time = Smooth(Accumulate(x,dt,nBinsX),smooth)';
	if pointProcess,
		map.z = map.count./(map.time+eps);
	else
		map.z = Smooth(Accumulate(x,z(:,2),nBinsX),smooth)';
		map.z = map.z./(map.count+eps);
	end
else
	% 2D (x and y)
	map.x = linspace(0,1,nBinsX);
	map.y = linspace(0,1,nBinsY);
	map.count = Smooth(Accumulate([x y],n,nBins),smooth)';
	map.time = Smooth(Accumulate([x y],dt,nBins),smooth)';
	if pointProcess,
		map.z = map.count./(map.time+eps);
	else
		map.z = Smooth(Accumulate([x y],z(:,2),nBins),smooth).';
		map.z = map.z./(map.count+eps);
	end
end

% Circular z
if strcmp(type,'circular'), map.z = wrap(angle(map.z),range); end

% Discard regions with insufficient sampling
map.z(map.time<=minTime) = 0;

