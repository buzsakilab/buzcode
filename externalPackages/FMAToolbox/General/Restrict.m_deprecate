function samples = Restrict(samples,intervals,varargin)

%Restrict - Keep only samples that fall in a given list of time intervals.
%
% Keep only samples (positions, spikes, LFP, etc.) that fall in a given list of
% time intervals.
%
% The remaining epochs can optionally be 'shifted' next to each other in time,
% removing the time gaps between them (which result from discarded samples),
% in which case they are also shifted globally to start at t = 0.
%
%  USAGE
%
%    samples = Restrict(samples,intervals,<options>)
%
%    samples        samples to restrict
%    intervals      list of (start,stop) pairs
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'shift'       shift remaining epochs together in time (default = 'off')
%    =========================================================================
%
%  NOTE
%
%    For more advanced time restriction of samples, use <a href="matlab:help InIntervals">InIntervals</a>.
%
%  SEE
%
%    See also ConsolidateIntervals, SubtractIntervals, ExcludeIntervals,
%    InIntervals, Restrict, FindInInterval, CountInIntervals, PlotIntervals.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
verbose = false;
shift = 'off';

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Restrict">Restrict</a>'' for details).');
end

if size(samples,1) == 1,
	samples = samples(:);
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Restrict">Restrict</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'shift',
			shift = varargin{i+1};
			if ~isstring_FMAT(shift,'on','off'),
				error('Incorrect value for property ''shift'' (type ''help <a href="matlab:help Restrict">Restrict</a>'' for details).');
			end

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Restrict">Restrict</a>'' for details).']);
	end
end

% Restrict
[status,interval,index] = InIntervals(samples,intervals);
samples = samples(status,:);

% Shift?
if strcmp(shift,'on'),
	% Discard interval IDs for samples which belong to none of the intervals
	interval = interval(status);
	% Samples in each interval will be shifted next to end of the previous interval
	% Let us call dt1 the time difference between interval 1 and interval 2. Interval 2 must be
	% shifted by dt1, interval 3 by dt1+dt2 (since interval 2 itself will be shifted by dt1),
	% interval 4 by dt1+dt2+dt3, etc.
	% 1) Compute the cumulative shifts dt1, dt1+dt2, dt1+dt2+dt3,...
	start = intervals(2:end,1); % beginning of next interval
	stop = intervals(1:end-1,2); % end of previous interval
	cumulativeShift = [0;cumsum(start-stop)];
	% 2) Assign these cumulative shifts to each sample
	shifts = cumulativeShift(interval);
	% 3) Shift
	samples(:,1) = samples(:,1) - shifts - intervals(1,1);
end
