function map = PhaseMap(positions,phases,varargin)

%PhaseMap - Compute spatial phase map.
%
%  Compute phase map (e.g. for a place cell), as well as occupancy and count maps.
%
%  USAGE
%
%    map = PhaseMap(positions,phases,<options>)
%
%    positions      position samples
%    phases         phases in radians
%    <options>      optional list of property-value pairs (see NOTE below)
%
%  NOTE
%
%    This function is provided for convenience. It simply calls <a href="matlab:help Map">Map</a>
%    using the same parameters. See this function for details about optional parameters.
%
%  OUTPUT
%
%    The outputs are the same as for <a href="matlab:help Map">Map</a>, except for map.z which is replaced by:
%
%    map.phase       average phase map (in radians)
%
%  SEE
%
%    See also Map, MapStats, FiringCurve, FiringMap, PlotColorMap.

% Copyright (C) 2009-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PhaseMap">PhaseMap</a>'' for details).');
end

% if size(samples,2) ~= 3,
%   error('Parameter ''positions'' is not a Nx3 matrix (type ''help <a href="matlab:help PhaseMap''">PhaseMap''</a> for details).');
% end

map = Map(positions,phases,varargin{:});
map.phase = map.z;
map = rmfield(map,'z');
