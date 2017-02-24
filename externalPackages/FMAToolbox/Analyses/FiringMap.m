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
%    <options>      optional list of property-value pairs (see NOTE below)
%
%  NOTE
%
%    This function is provided for convenience. It simply calls <a href="matlab:help Map">Map</a> and
%    <a href="matlab:help MapStats">MapStats</a> using the same parameters. See these functions for details
%    about optional parameters.
%
%  OUTPUT
%
%    The outputs are the same as for <a href="matlab:help Map">Map</a> and <a href="matlab:help MapStats">MapStats</a>, except for map.z which
%    is replaced by:
%
%    map.rate       average firing rate map (in Hz)
%
%  SEE
%
%    See also Map, MapStats, FiringCurve, PlotColorMap.

% Copyright (C) 2004-2011 by Michaël Zugaro
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
  warning('Parameter ''positions'' is not a Nx3 matrix - using the red light.');
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
		case {'center','threshold','minsize','type'},
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
