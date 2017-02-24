%FindField - Find field in a map, e.g. firing field in a firing map.
%
% Find field in a map, i.e. the connex area around the peak, where the values
% are above a given threshold.
%
% Note: This function should not be called directly. Use <a href="matlab:help FiringMap">FiringMap</a>, <a href="matlab:help FiringCurve">FiringCurve</a>
% or <a href="matlab:help MapStats">MapStats</a> instead.
%
%  USAGE
%
%    field = FindField(map,x,y,threshold,circX,circY)
%
%    map            firing map
%    x              abscissa (column) of the peak
%    y              ordinate (row) of the peak
%    threshold      min firing rate within the field
%    circX          true if x coordinate is circular
%    circY          true if y coordinate is circular
%
%  SEE
%
%    See also DefineField, MapStats, FiringMap, FiringCurve, IsInZone,
%    PlotColorMap.

% Copyright (C) 2004-2012 by Michaël Zugaro
