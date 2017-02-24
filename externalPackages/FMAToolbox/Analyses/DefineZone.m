function zone = DefineZone(s,shape,points)

%DefineZone - Define a restricted zone to which analyses can be circumscribed.
%
% Define a restricted zone in the environment, to which analyses can subsequently be
% circumscribed using <a href="matlab:help IsInZone">IsInZone</a>. Multiple zones can easily be combined (see example).
%
%  USAGE
%
%    zone = DefineZone(size,shape,points)
%
%    size           environment size [h v]
%    shape          'rectangle' or 'circle'
%    points         list of points defining the field
%                     rectangle: [x y width height]
%                     circle: [x y radius]
%
%  EXAMPLE
%
%    % Define a rectangular and a circular zones, and combine them
%    r = DefineZone([100 100],'rectangle',[20 30 50 20]);
%    c = DefineZone([100 100],'circle',[50 70 15]);
%    combined = c|r;
%
%  SEE
%
%    See also FiringMap, IsInZone.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 3,
  error('Incorrect number of parameters (type ''help <a href="matlab:help DefineZone">DefineZone</a>'' for details).');
end
if ~isdvector(s,'#2','>0'),
  error('Incorrect size (type ''help <a href="matlab:help DefineZone">DefineZone</a>'' for details).');
end
if ~isstring_FMAT(shape,'rectangle','circle'),
  error('Incorrect shape (type ''help <a href="matlab:help DefineZone">DefineZone</a>'' for details).');
end
points = round(points);

width = s(1);
height = s(2);

zone = logical(zeros(width,height));

if strcmp(lower(shape),'rectangle'),
	if length(points) ~= 4,
		error('Incorrect points (type ''help <a href="matlab:help DefineZone">DefineZone</a>'' for details).');
	end
	x = points(1);
	y = points(2);
	w = points(3);
	h = points(4);
	zone(y:y+h-1,x:x+w-1) = 1;
elseif strcmp(lower(shape),'circle'),
	if length(points) ~= 3,
		error('Incorrect points (type ''help <a href="matlab:help DefineZone">DefineZone</a>'' for details).');
	end
	x = points(1);
	y = points(2);
	r = points(3);
	for i = 1:width,
		for j = 1:height,
			if (x-i)^2+(y-j)^2 <= r^2,
				zone(j,i) = 1;
			end
		end
	end
end

