function [up,down] = ZeroCrossings(data,varargin)

%ZeroCrossings - Test zero crossings in a given time series.
%
%  This assumes a minimum of 10 samples per positive/negative phase.
%
%  USAGE
%
%    [up,down] = ZeroCrossings(samples,<options>)
%
%    samples        an Nx2 matrix of (timestamp,value) pairs
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      standard deviation for Gaussian kernel (default = 0)
%    =========================================================================
%
%  OUTPUT
%
%    up             logical indices indicating upward zero crossings
%    down           logical indices indicating downward zero crossings
%
%  NOTE
%
%    This finds the points in the signal closest to the zero crossings
%    (the actual zero crossings may not appear in the time series).
%    Hence, if the signal changes very abruptly, e.g. if X(t)=A
%    and X(t+1)=-A (where A is some large constant), this function
%    will consider index t as a downward zero crossing, although
%    the value A is very different from 0.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ZeroCrossings">ZeroCrossings</a>'' for details).');
end

if size(data,2) ~= 2,
  error('Parameter ''data'' is not a Nx2 matrix (type ''help <a href="matlab:help ZeroCrossings">ZeroCrossings</a>'' for details).');
end

% Defaults
smooth = 0;

% Parse parameter list
for j = 1:2:length(varargin),
	if ~ischar(varargin{j}),
		error(['Parameter ' num2str(j+7) ' is not a property (type ''help <a href="matlab:help ZeroCrossings">ZeroCrossings</a>'' for details).']);
	end
	switch(lower(varargin{j})),
		case 'smooth',
			smooth = varargin{j+1};
			if ~isdvector(smooth,'>=0') || length(smooth) > 2,
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help ZeroCrossings">ZeroCrossings</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{j}) ''' (type ''help <a href="matlab:help ZeroCrossings">ZeroCrossings</a>'' for details).']);
	end
end

data(:,2) = Smooth(data(:,2),smooth);

% Find downward and upward going zero-crossings
previous = data(1:end-1,2);
current = data(2:end,2);
down = previous > 0 & current < 0;down(end+1) = 0;
up = previous < 0 & current > 0;up(end+1) = 0;

