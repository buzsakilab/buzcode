function dt = DistanceTransform(b)

%DistanceTransform - Compute distance transform of a binary matrix (distance to nearest 1 value).
%
%  This function assigns to every point (x,y) in a binary matrix B the distance
%  to the nearest point in B with value 1. It uses the Euclidean metric for
%  computing distances. The code is based on the algorithm by Meijster et al.
%  (2000).
%
%  USAGE
%
%    dt = DistanceTransform(b)
%
%    b              binary matrix
%

% Copyright (C) 2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% The code is not self explanatory, but is a direct transcription of the algorithms provided
% in Meijster et al. (2000)

[m,n] = size(b);
g = nan(m,n);

for x = 1:n,
	% Scan 1
	if b(1,x),
		g(1,x) = 0;
	else
		g(1,x) = m + n;
	end
	for y = 2:m,
		if b(y,x),
			g(y,x) = 0;
		else
			g(y,x) = 1 + g(y-1,x);
		end
	end
	% Scan 2
	for y = m-1:-1:1,
		if g(y+1,x) < g(y,x),
			g(y,x) = 1 + g(y+1,x);
		end
	end
end

for y = 1:n,
	q = 1;
	s(1) = 1;
	t(1) = 1;
	for u = 2:n,
		% Scan 3
		while true,
			x = t(q);
			i = s(q);
			f1 = (x-i)^2 + g(y,i)^2;
			i = u;
			f2 = (x-i)^2 + g(y,i)^2;
			if f1 <= f2, break; end
			q = q - 1;
			if q < 1, break; end
		end
		if q < 1,
			q = 1;
			s(1) = u;
		else
			i = s(q);
			a = (u^2-i^2+g(y,u)^2-g(y,i)^2);
			b = 2*(u-i);
			w = 1 + floor(a/b);
			if w <= n,
				q = q + 1;
				s(q) = u;
				t(q) = w;
			end
		end
	end
	for u = n:-1:1,
		% Scan 4
		x = u;
		i = s(q);
		dt(y,u) = sqrt((x-i)^2 + g(y,i)^2);
		if u == t(q),
			q = q - 1;
		end
	end
end
