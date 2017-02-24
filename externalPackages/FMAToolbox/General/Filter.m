function filtered = Filter(samples,varargin)

%Filter - Filter samples.
%
%  USAGE
%
%    filtered = Filter(samples,<options>)
%
%    samples        samples given as a list of (t,v1,v2,v3...) tuples
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'passband'    pass frequency range
%     'stopband'    stop frequency range
%     'order'       filter order (default = 4)
%     'ripple'      filter ripple (default = 20)
%     'nyquist'     nyquist frequency (default = 625)
%     'filter'      choose filter type between 'cheby2' (default) and 'fir1'
%    =========================================================================
%

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
passband = [];
stopband = [];
order = 4;
ripple = 20;
nyquist = 625;
type = 'cheby2';

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
end

% Check parameter sizes
if size(samples,2) < 2,
	error('Parameter ''samples'' is not a matrix (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help Filter">Filter</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'passband',
			if ~isempty(stopband),
				error('Cannot specify both a passband and stopband (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end
			passband = varargin{i+1};
			if ~isdvector(passband,'#2','>=0'),
				error('Incorrect value for ''passband'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end

		case 'stopband',
			if ~isempty(passband),
				error('Cannot specify both a passband and stopband (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end
			stopband = varargin{i+1};
			if ~isdvector(stopband,'#2','>=0'),
				error('Incorrect value for ''stopband'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end

		case 'filter',
			type = lower(varargin{i+1});
			if ~isstring_FMAT(type,'cheby2','fir1'),
				error(['Unknown filter type ''' type ''' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).']);
			end

		case 'order',
			order = lower(varargin{i+1});
			if ~isiscalar(order,'>0'),
				error('Incorrect value for ''order'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end

		case 'ripple',
			ripple = lower(varargin{i+1});
			if ~isiscalar(ripple,'>0'),
				error('Incorrect value for ''ripple'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end

		case 'nyquist',
			nyquist = varargin{i+1};
			if ~isiscalar(nyquist,'>0'),
				error('Incorrect value for property ''nyquist'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).']);

	end
end

if isempty(passband) && isempty(stopband),
	error('Missing passband or stopband (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
end

switch(type),
	case 'cheby2',
		if ~isempty(passband),
			if passband(1) == 0,
				[b a] = cheby2(order,ripple,passband(2)/nyquist,'low');
			else
				[b a] = cheby2(order,ripple,passband/nyquist);
			end
		else
			[b a] = cheby2(order,ripple,stopband/nyquist,'stop');
		end
	case 'fir1',
		if ~isempty(passband),
			if passband(1) == 0,
				[b a] = fir1(order,passband(2)/nyquist,'low');
			else
				[b a] = fir1(order,passband/nyquist);
			end
		else
			[b a] = fir1(order,stopband/nyquist,'stop');
		end
end
filtered(:,1) = samples(:,1);
for i = 2:size(samples,2),
	filtered(:,i) = filtfilt(b,a,samples(:,i));
end
