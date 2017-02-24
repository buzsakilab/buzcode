function [periods,in] = Threshold(x,criterion,threshold,varargin)

%Threshold - Find periods above/below threshold.
%
%  Find periods where samples lay above or below a given threshold. Optionally,
%  keep only periods of sufficient duration and ignore brief interruptions.
%
%  USAGE
%
%    [periods,in] = Threshold(x,threshold,<options>)
%
%    x              list of time-value pairs
%    criterion      one of '>', '>=', '<' or '<='
%    threshold      threshold
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'min'         minimum duration for inclusion (default = 0)
%     'max'         maximum duration of ignored interruptions (default = 0)
%    =========================================================================
%
%  OUTPUT
%
%    periods        list of [start stop] pairs
%    in             list of [t s] pairs, where s is 1 if the the values
%                   match the inclusion criterion (and 0 otherwise)
%
%  SEE
%
%    See also QuietPeriods.

% Copyright (C) 2008-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
m = 0;
M = 0;

if nargin < 3,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Threshold">Threshold</a>'' for details).');
end
if ~isstring_FMAT(criterion,'>','>=','<','<='),
	error('Incorrect criterion (type ''help <a href="matlab:help Threshold">Threshold</a>'' for details).');
end

if mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Threshold">Threshold</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help Threshold">Threshold</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'min',
			m = varargin{i+1};
			if ~isdscalar(m,'>=0'),
				error('Incorrect value for property ''min'' (type ''help <a href="matlab:help Threshold">Threshold</a>'' for details).');
			end
		case 'max',
			M = varargin{i+1};
			if ~isdscalar(m,'>=0'),
				error('Incorrect value for property ''min'' (type ''help <a href="matlab:help Threshold">Threshold</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Threshold">Threshold</a>'' for details).']);
	end
end

% Determine beginning/end of included periods
ok = eval(['x(:,2)' criterion 'threshold']);
crossings = diff(ok); % yields -1 for in->out crossings, and 1 for out->in crossings
start = find(crossings == 1);
stop = find(crossings == -1);

% The previous code would ignore periods beginning at the first sample, or ending at the last sample; correct for this
if ok(1),
	start = [1;start];
end
if ok(end),
	stop = [stop;length(ok)];
end

% Determine durations of excluded periods, and include brief ones
durations = x(start(2:end),1) - x(stop(1:end-1),1);
ignore = find(durations <= M);
start(ignore+1) = [];
stop(ignore) = [];

% Keep only long enough periods
durations = x(stop,1)-x(start,1);
discard = durations < m;
start(discard) = [];
stop(discard) = [];

% Outputs
periods = [x(start,1) x(stop,1)];
in = logical(zeros(size(x,1),1));
for i = 1:length(start),
	in(start(i):stop(i)) = 1;
end
