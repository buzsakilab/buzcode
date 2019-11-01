function [periods,movement] = MovementPeriods(v,velocity,duration,brief)

%QuietPeriods - Find periods of immobility.
%
%  Find periods of movement, i.e. periods of sufficient duration
%  where instantaneous linear velocity is above threshold. Brief movements
%  can be ignored.
%
%  USAGE
%
%    [periods,movement] = MovementPeriods(v,velocity,duration,brief)
%
%    v              linear velocity samples [t v]
%    velocity       minimun velocity of a movement period
%    duration       minimum duration of a movement period
%    brief          optional maximum duration of a 'brief' movement
%
%  OUTPUT
%
%    periods        list of [start stop] pairs
%    movement     list of [t s] pairs, where s is 1 if the animal
%                   is quiet at time t (and 0 otherwise)
%
%
%    Modified from QuietPeriods.m
%    AntonioFR, 03/2019


if nargin < 3,
	error('Incorrect number of parameters (type ''help <a href="matlab:help QuietPeriods">QuietPeriods</a>'' for details).');
end
if nargin == 3,
	brief = 0;
end

% Determine beginning/end of mov periods
above = v(:,2) > velocity;
crossings = diff(above); % yields -1 for upward crossings, and 1 for downward crossings
start = find(crossings == 1);
stop = find(crossings == -1);

% The previous code would ignore mov periods beginning at the first sample, or ending at the last sample; correct for this
if above(1),
	start = [1;start];
end
if above(end),
	stop = [stop;length(above)];
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
movement = zeros(size(v,1),2);
movement(:,1) = v(:,1);
for i = 1:length(start),
	movement(start(i):stop(i),2) = 1;
end
