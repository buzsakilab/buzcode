function distance = Distance(positions,reference,varargin)

%Distance - Compute instantaneous distance to a reference point.
%
%  USAGE
%
%    distance = Distance(positions,reference,<options>)
%
%    positions      a list of position samples
%    reference      reference point (same coordinates as positions)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'type'        'l' if X is linear (default), 'c' if X is circular (for
%                    1D positions only, which should be in [0..1])
%    =========================================================================
%
%  SEE
%
%    See also Diff.

% Copyright (C) 2004-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
type = 'linear';

if nargin < 2 | mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Distance">Distance</a>'' for details).');
end
if ~issamples(positions,'#1') && ~issamples(positions,'#2'),
	error('Incorrect positions - should be a vector (type ''help <a href="matlab:help Distance">Distance</a>'' for details).');
end
if ~isdvector(reference,'#1') && ~isdvector(reference,'#2'),
	error('Incorrect reference - should be a vector (type ''help <a href="matlab:help Distance">Distance</a>'' for details).');
end
if length(reference) ~= size(positions,2)-1,
	error('Positions and reference have incompatible sizes (type ''help <a href="matlab:help Distance">Distance</a>'' for details).');
end

% Parse parameter list
for j = 1:2:length(varargin),
	if ~ischar(varargin{j}),
		error(['Parameter ' num2str(j+2) ' is not a property (type ''help <a href="matlab:help Distance">Distance</a>'' for details).']);
	end
	switch(lower(varargin{j})),
		case 'type',
			type = varargin{j+1};
			if ~isstring_FMAT(type,'l','c'),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help Distance">Distance</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{j}) ''' (type ''help <a href="matlab:help Distance">Distance</a>'' for details).']);
	end
end

distance = [];
if isempty(positions), return; end

distance = positions(:,1);

if issamples(positions,'#1'),
	if type(1) == 'l',
		distance(:,2) = abs(positions(:,2)-reference(1));
	else
		% Make sure X is normalized
		if max(positions(:,2)) > 1 || min(positions(:,2)) < 0,
			positions(:,2) = ZeroToOne(positions(:,2));
			warning('Positions should contain values in [0 1]. The data will now be transformed accordingly.');
		end
		distance1 = abs(positions(:,2)-reference(1));
		distance2 = 1-abs(positions(:,2)-reference(1));
		distance(:,2) = min([distance1 distance2],[],2);
	end
else
	distance(:,2) = sqrt((positions(:,2)-reference(1)).^2+(positions(:,3)-reference(2)).^2);
end