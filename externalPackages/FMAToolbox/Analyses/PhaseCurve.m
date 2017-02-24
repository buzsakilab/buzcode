function map = PhaseCurve(samples,phases,varargin)

%PhaseCurve - Compute spatial phase curve.
%
%  Compute phase curve, as well as occupancy and count curves.
%
%  USAGE
%
%    curve = PhaseCurve(samples,phases,<options>)
%
%    samples        e.g. angles or linear positions across time
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
%    The outputs are the same as for <a href="matlab:help Map">Map</a>, except for curve.z which is replaced by:
%
%    curve.phase       average phase curve (in radians)
%
%  SEE
%
%    See also PhaseMap, Map, MapStats, FiringCurve, FiringMap, PlotXY.

% Copyright (C) 2009-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PhaseCurve">PhaseCurve</a>'' for details).');
end

if ~isdmatrix(samples) | size(samples,2) ~= 2,
  error('Parameter ''samples'' is not a Nx2 matrix (type ''help <a href="matlab:help PhaseCurve">PhaseCurve</a>'' for details).');
end

map = Map(samples,phases,'type','circular',varargin{:});
map.phase = map.z;
map = rmfield(map,'z');
