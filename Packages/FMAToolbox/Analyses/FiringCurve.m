function [curve,stats] = FiringCurve(samples,spikes,varargin)

%FiringCurve - Compute firing curve (e.g. for a head direction cell).
%
%  Compute firing curve (e.g. for a head direction cell, or a place or grid cell
%  on a linear track), as well as occupancy and spike count curves. Optionnally,
%  curve statistics can also be computed, including in-field peak and mean firing
%  rates, firing field width, etc.
%
%  USAGE
%
%    [curve,stats] = FiringCurve(samples,spikes,<options>)
%
%    samples        e.g. angles or linear samples across time
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
%    curve.rate      average firing rate curve (in Hz)
%
%  SEE
%
%    See also Map, MapStats, FiringMap, PlotXY.

% Copyright (C) 2005-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help FiringCurve">FiringCurve</a>'' for details).');
end

if size(samples,2) ~= 2,
  error('Parameter ''samples'' is not a Nx2 matrix (type ''help <a href="matlab:help FiringCurve">FiringCurve</a>'' for details).');
end

im = 1;argsm = {};
is = 1;argss = {};
% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FiringCurve">FiringCurve</a>'' for details).']);
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

curve = Map(samples,spikes,argsm{:});
if nargout == 2,
	stats = MapStats(curve,argss{:});
end
curve.rate = curve.z;
curve = rmfield(curve,'z');
