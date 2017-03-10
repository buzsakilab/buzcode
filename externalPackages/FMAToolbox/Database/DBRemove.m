function DBRemove(database)

%DBRemove - Remove database.
%
%  Actual removal will require confirmation from the user.
%
%  USAGE
%
%    DBRemove(database)
%
%    database       database name
%
%  SEE
%
%    See also DBRemove, DBRemoveFigures, DBRemoveVariables.
%

% Copyright (C) 2007-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

% Check parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help DBRemove">DBRemove</a>'' for details).');
end
if ~isstring(database),
  error('Incorrect database name (type ''help <a href="matlab:help DBRemove">DBRemove</a>'' for details).');
end

% Select database
try
	h = mym(['use ' database]);
catch
	error(['Could not find database ''' database '''.']);
end

% Print warning and display # figures and # variables to remove
disp(' ');
f = mym(['select name from figures']);
if length(f) > 1, sf = 's'; else sf = ''; end
v = mym(['select name from variables']);
if length(v) > 1, sv = 's'; else sv = ''; end
disp(['This would remove ''' database ''' (' int2str(length(f.name)) ' figure' sf ', ' int2str(length(v.name)) ' variable' sv ').']);

% Confirm
s = lower(input('Type ''remove'' to confirm: ','s'));
if ~strcmp(s,'remove'),
	disp('*** Cancelled ***');
	return
end

% Remove
mym(['drop database ' database]);
disp('Database removed.');
