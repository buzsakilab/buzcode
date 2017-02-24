function [x,m,e] = RunningAverage(X,Y,varargin)

%RunningAverage - Compute running linear or angular average.
%
% Computes the running average of y=f(x). Variable y can be linear or circular
% (use radians). The error bars are standard errors of the mean for linear
% data, or 95% confidence intervals for circular data.
%
%  USAGE
%
%    [x,m,e] = RunningAverage(x,y,<options>)
%
%    x              x variable (e.g., time)
%    y              values at x
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'window'      averaging window size (default (max-min)/10)
%     'overlap'     overlap between successive windows (default = 0.8*window)
%     'limits'      x limits to use instead of min and max
%     'type'        either 'linear' or 'circular' (default 'linear')
%    =========================================================================
%
%  OUTPUT
%
%    x              new x variable
%    m              running average
%    e              sem for linear variables, otherwise 95% confidence intervals

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
type = 'linear';
limits = [min(X) max(X)];
nBins = 10;
window = 0;
overlap = [];

if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help RunningAverage">RunningAverage</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help RunningAverage">RunningAverage</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'type',
			type = varargin{i+1};
			if ~isstring_FMAT(type,'linear','circular'),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help RunningAverage">RunningAverage</a>'' for details).');
			end

		case 'window',
			window = varargin{i+1};
			if ~isdscalar(window,'>0'),
				error('Incorrect value for property ''window'' (type ''help <a href="matlab:help RunningAverage">RunningAverage</a>'' for details).');
			end

		case 'overlap',
			overlap = varargin{i+1};
			if ~isdscalar(overlap,'>=0'),
				error('Incorrect value for property ''overlap'' (type ''help <a href="matlab:help RunningAverage">RunningAverage</a>'' for details).');
			end

		case 'limits',
			limits = varargin{i+1};
			if ~isdvector(limits,'#2','<'),
				error('Incorrect value for property ''limits'' (type ''help <a href="matlab:help RunningAverage">RunningAverage</a>'' for details).');
			end

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help RunningAverage">RunningAverage</a>'' for details).']);

	end
end

x0 = limits(1);
x1 = limits(2);
if window == 0,
	window = (x1-x0)/nBins;
end
if isempty(overlap),
	overlap = 0.8*window;
end

% Loop through data
i = 1;
xi = x0+window/2;
while xi+window/2 <= x1,
	x(i,1) = xi;
	ok = InIntervals(X,xi+[-0.5 0.5]*window);
	if sum(ok) == 0,
		m(i,1) = NaN;
		e(i,:) = [NaN NaN];
	elseif strcmp(type,'circular'),
		[M,C] = CircularConfidenceIntervals(Y(ok));
		m(i,1) = M;
		e(i,:) = C;
	else
		m(i,1) = nanmean(Y(ok));
		n = sum(ok);
		s = nanstd(Y(ok))/sqrt(n);
		e(i,1) = m(i)-s;
		e(i,2) = m(i)+s;
	end
	xi = xi + window - overlap;
	i = i + 1;
end
