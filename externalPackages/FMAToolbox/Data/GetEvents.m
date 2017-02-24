function events = GetEvents(selection,varargin)

%GetEvents - Get events.
%
%  USAGE
%
%    events = GetEvents(selection,<options>)
%
%    selection      optional list of event descriptions; see examples below
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'output'      'times' returns the event times, 'indices' returns a
%                   list of matching row indices, 'logical' returns a
%                   list of matching row logical indices, and 'descriptions'
%                   returns a list of descriptions (default = 'times')
%    =========================================================================
%
%  EXAMPLES
%
%    % Get all event descriptions
%    evt = GetEvents('output','descriptions');
%
%    % Get all event times
%    evt = GetEvents;
%
%    % Get timestamps of events with description 'Ripple' or 'Sharp Wave'
%    evt = GetEvents({'Ripple','Sharp Wave'});
%
%    % Get indices of events with description 'Ripple'
%    idx = GetEvents({'Ripple'},'output','indices');
%    % When there is only one regexp, the {} can be omitted:
%    idx = GetEvents('Ripple','output','indices');
%
%    % Get logical indices of events with description starting with any letter
%    % but 'R' and ending with ' beginning'
%    m = GetEvents({'[^R].* beginning'},'output','logical');
%
%  NOTE
%
%    Type 'help regexp' for details on pattern matching using regular expressions.

% Copyright (C) 2004-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA;
if isempty(DATA),
	error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end

% Default values
output = 'times';

if nargin == 0,
	selection = '.*';
end

if mod(length(varargin),2) ~= 0,
	varargin = {selection varargin{:}};
	selection = [];
else
	if isa(selection,'char'), selection = {selection}; end
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+firstIndex) ' is not a property (type ''help <a href="matlab:help GetEvents">GetEvents</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'output',
			output = lower(varargin{i+1});
			if ~isstring_FMAT(output,'times','indices','logical','descriptions'),
				error('Incorrect value for property ''output'' (type ''help <a href="matlab:help GetEvents">GetEvents</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help GetEvents">GetEvents</a>'' for details).']);
	end
end

events = DATA.events;
if isempty(events.time), return; end

% Selected events only
nPatterns = length(selection);
if nPatterns == 0,
	selected = logical(ones(size(events.time)));
else
	selected = logical(zeros(size(events.time)));
	for i = 1:nPatterns,
		pattern = ['^' selection{i} '$'];
		matches = GetMatches(regexp(events.description,pattern));
		selected = selected | matches;
	end
end
if strcmp(output,'times'),
	events = events.time(selected,:);
elseif strcmp(output,'indices'),
	events = find(selected);
elseif strcmp(output,'logical'),
	events = selected;
elseif strcmp(output,'descriptions'),
	events = unique({events.description{selected}})';
end

function m = GetMatches(c)

m = zeros(length(c),1);
for i = 1:length(c),
	m(i) = ~isempty(c{i});
end
