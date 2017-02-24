function [map,stats] = FiringMap(positions,spikes,varargin)

%FiringMap - Compute firing map (e.g. for a place cell).
%
%  Compute firing map (e.g. for a place cell), as well as occupancy and spike
%  count maps. Optionnally, field statistics can also be computed, including
%  in-field peak and mean firing rates, firing field size, etc.
%
%  USAGE
%
%    [map,stats] = FiringMap(positions,spikes,<options>)
%
%    positions      position samples
%    spikes         spike timestamps
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
%     'type'        two letters (one for X and one for Y) indicating which
%                   coordinates are linear ('l') and which are circular ('c')
%                   (default 'll')
%    =========================================================================
%
%  OUTPUT
%
%    map.x               x bins
%    map.y               y bins
%    map.rate            average firing rate map (in Hz)
%    map.count           spike count map
%    map.time            occupancy map (in s)
%
%    stats.x             abscissa of the peak rate (in bins)
%    stats.y             ordinate of the peak rate (in bins)
%    stats.peak          in-field peak rate
%    stats.mean          in-field mean value
%    stats.size          field size (in bins)
%    stats.field         field (1 = bin in field, 0 = bin not in field)
%    stats.fieldX        field x boundaries (in bins)
%    stats.fieldY        field y boundaries (in bins)
%    stats.specificity   spatial specificity (Skaggs et al., 1993)
%
%  NOTE
%
%    This function is provided for convenience. It simply calls <a href="matlab:help Map">Map</a> and <a href="matlab:help MapStats">MapStats</a>
%    using the same parameters. The outputs are the same except for map.z which
%    is replaced by map.rate.

%  SEE
%
%    See also Map, MapStats, FiringCurve, PlotColorMap.

% Copyright (C) 2004-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help FiringMap">FiringMap</a>'' for details).');
end

if size(positions,2) ~= 3,
  warning('Parameter ''positions'' is not a Nx3 matrix - using the first light.');
  positions = positions(:,1:3);
end

im = 1;argsm = {};
is = 1;argss = {};
% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FiringMap">FiringMap</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'type',
			argss{is} = 'type';
			argss{is+1} = varargin{i+1}(1:2);
			is = is+2;
			argsm{im} = 'type';
			argsm{im+1} = [varargin{i+1}(1:2) 'l'];
			im = im+2;
		case {'threshold','minsize','minpeak','debug'},
			argss{is} = varargin{i};
			argss{is+1} = varargin{i+1};
			is = is+2;
		otherwise,
			argsm{im} = varargin{i};
			argsm{im+1} = varargin{i+1};
			im = im+2;
  end
end

map = Map(positions,spikes,argsm{:});
if nargout == 2,
	stats = MapStats(map,argss{:});
end
map.rate = map.z;
map = rmfield(map,'z');
