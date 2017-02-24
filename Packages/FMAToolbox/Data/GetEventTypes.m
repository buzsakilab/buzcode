function descriptions = GetEventTypes(selection)

%GetEventTypes - Get event types (= unique descriptions).
%
%  USAGE
%
%    descriptions = GetEventTypes(selection)
%
%    selection      optional patterns of event descriptions; see examples below
%
%  EXAMPLES
%
%    % Get all event descriptions
%    descriptions = GetEventTypes();
%
%    % Get events with description starting with any letter
%    % but 'R' and ending with ' beginning'
%    m = GetEventTypes({'[^R].* beginning'});
%
%  NOTE
%
%    Type 'help regexp' for details on pattern matching using regular expressions.
%
%  SEE
%
%    See also GetEvents.


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

if nargin == 0,
	selection = '.*';
end
if isa(selection,'char'), selection = {selection}; end

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
descriptions = unique({events.description{selected}})';

function m = GetMatches(c)

m = zeros(length(c),1);
for i = 1:length(c),
	m(i) = ~isempty(c{i});
end
