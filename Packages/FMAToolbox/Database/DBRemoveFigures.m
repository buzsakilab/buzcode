function DBRemoveFigures(query)

%DBRemoveFigures - Remove all figures that match given criteria.
%
%  Actual removal will require confirmation from the user.
%
%  USAGE
%
%    DBRemoveFigures(query)
%
%    query          optional figure list query (WHERE clause; see Example)
%
%  SEE
%
%    See also DBAddFigure, DBGetFigures.
%

% Copyright (C) 2007-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

% Optional query provided?
if nargin == 0,
	query = '';
end

% Edit query
query = strtrim(query);
query = regexprep(query,'^where','');
if ~isempty(query), query = [' where ' query]; end

% Get database name
database = DBUse;

% Query database
f = mym(['select eid,name from figures' query]);

% Make sure query results are not empty
if isempty(f.eid),
	if isempty(query),
		warning(['No figures in ''' database '''.']);
	else
		warning(['No figures match (' query ').']);
	end
	return
end

% Print warning and display figures to remove
disp(' ');
disp(['This would remove the following figures from ''' database ''':']);
DBDisplay(f);

% Confirm
s = lower(input('Type ''remove'' to confirm: ','s'));
if ~strcmp(s,'remove'),
	disp('*** Cancelled ***');
	return
end

% Remove
mym(['delete from figures' query]);
disp('Figures removed.');
