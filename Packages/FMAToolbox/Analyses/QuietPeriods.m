function [periods,quiescence] = QuietPeriods(v,velocity,duration,brief)

%QuietPeriods - Find periods of immobility.
%
%  Find periods of immobility, i.e. periods of sufficient duration
%  where instantaneous linear velocity remains low. Brief movements
%  can be ignored.
%
%  USAGE
%
%    [periods,quiescence] = QuietPeriods(v,velocity,duration,brief)
%
%    v              linear velocity samples [t v]
%    velocity       maximum velocity of a quiet period
%    duration       minimum duration of a quiet period
%    brief          optional maximum duration of a 'brief' movement
%
%  OUTPUT
%
%    periods        list of [start stop] pairs
%    quiescence     list of [t s] pairs, where s is 1 if the animal
%                   is quiet at time t (and 0 otherwise)
%
%  SEE
%
%    See also BrainStates, PlotIntervals.

% Copyright (C) 2008-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 3,
	error('Incorrect number of parameters (type ''help <a href="matlab:help QuietPeriods">QuietPeriods</a>'' for details).');
end
if nargin == 3,
	brief = 0;
end

% Determine beginning/end of quiet periods
below = v(:,2) < velocity;
crossings = diff(below); % yields -1 for upward crossings, and 1 for downward crossings
start = find(crossings == 1);
stop = find(crossings == -1);

% The previous code would ignore quiet periods beginning at the first sample, or ending at the last sample; correct for this
if below(1),
	start = [1;start];
end
if below(end),
	stop = [stop;length(below)];
end

% Determine durations of movements, and discard brief ones
durations = v(start(2:end),1) - v(stop(1:end-1),1);
ignore = find(durations <= brief);
start(ignore+1) = [];
stop(ignore) = [];

% Keep only long enough periods
durations = v(stop,1)-v(start,1);
discard = durations < duration;
start(discard) = [];
stop(discard) = [];

% Outputs
periods = [v(start,1) v(stop,1)];
quiescence = zeros(size(v,1),2);
quiescence(:,1) = v(:,1);
for i = 1:length(start),
	quiescence(start(i):stop(i),2) = 1;
end
