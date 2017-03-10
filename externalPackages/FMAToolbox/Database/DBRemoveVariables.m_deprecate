function DBRemoveVariables(query)

%DBRemoveVariables - Remove all variables that match given criteria.
%
%  Actual removal will require confirmation from the user.
%
%  USAGE
%
%    DBRemoveVariables(query)
%
%    query          optional figure list query (WHERE clause; see Example)
%
%  SEE
%
%    See also DBAddVariable, DBGetVariables.
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
f = mym(['select eid,name from variables' query]);

% Make sure query results are not empty
if isempty(f.eid),
	if isempty(query),
		warning(['No variables in ''' database '''.']);
	else
		warning(['No variables match (' query ').']);
	end
	return
end

% Print warning and display variables to remove
disp(' ');
disp(['This would remove the following variables from ''' database ''':']);
DBDisplay(f);

% Confirm
s = lower(input('Type ''remove'' to confirm: ','s'));
if ~strcmp(s,'remove'),
	disp('*** Cancelled ***');
	return
end

% Remove
mym(['delete from variables' query]);
disp('Variables removed.');
