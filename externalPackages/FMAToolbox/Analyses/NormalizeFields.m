function normalized = NormalizeFields(fields,varargin)

%NormalizeFields - Normalize one or more firing fields in space and rate.
%
%  For each field, interpolate space in [0..1] and normalize z values
%  (typically, firing rates for place fields). Z values outside the fields
%  should be set to zero (see example below).
%
%  Current implementation is only for 1D environments.
%
%  USAGE
%
%    normalized = NormalizeFields(fields,<options>)
%
%    fields         firing fields (MxN: M fields, N bins)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'rate'        'off' to disable rate normalization (default = 'on')
%    =========================================================================
%
%  EXAMPLE
%
%    [c1,s1] = FiringCurve(positions,spikes1);
%    c1.rate(~c1.field) = 0;
%    c2 = FiringCurve(positions,spikes2);
%    c2.rate(~c2.field) = 0;
%    n = NormalizeFields([c1.rate;c2.rate]);
%    average = mean(n);
%    error = std(n);

% Copyright (C) 2012-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
rate = 'on';

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help NormalizeFields">NormalizeFields</a>'' for details).');
end

% Check parameter sizes
if ~isdmatrix(fields),
	error('Firing fields should be MxN matrices (type ''help <a href="matlab:help NormalizeFields">NormalizeFields</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help NormalizeFields">NormalizeFields</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'rate',
			rate = varargin{i+1};
			if ~isstring_FMAT(rate,'on','off'),
				error('Incorrect value for property ''rate'' (type ''help <a href="matlab:help NormalizeFields">NormalizeFields</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help NormalizeFields">NormalizeFields</a>'' for details).']);
	end
end

[m,n] = size(fields);

% First, align fields to first bin

% For each field, find start index, i.e. transition from 0 to 1 (0: outside field, 1: inside)
f = logical(fields);
transitions = [zeros(m,1) diff(f,1,2)];
[i,j] = find(transitions==1);
start = ones(m,1);
start(i) = j;
% Shift to align left
aligned = CircularShift(fields,-(start-1));

% Second, 'spread' fields across all bins

% For each field, find stop index, i.e. transition from 1 to 0
f = logical(aligned);
transitions = [zeros(m,1) diff(f,1,2)];
[i,j] = find(transitions==-1);
stop = n*ones(m,1);
stop(i) = j-1;
% Interpolate
normalized = zeros(m,n);
for i = 1:m,
	normalized(i,:) = interp1(aligned(i,1:stop(i)),linspace(1,stop(i),n));
end

% Third, normalize z values
if strcmp(rate,'on'),
	M = max(normalized,[],2);
	normalized = normalized ./ repmat(M,1,n);
end
