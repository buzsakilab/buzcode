function database = DBUse(database)

%DBUse - Set (or determine) current database.
%
% Set (or determine) current database to store/retrieve processed data.
%
%  USAGE
%
%    DBUse(database)
%    database = DBUse()
%
%    database           database name
%
%  SEE
%
%    See also DBConnect.
%

% Copyright (C) 2007-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

if nargin == 0,
	try
		h = mym('select database();');
	catch
		error(['Could not determine current database (check DB server connection).']);
	end
	database = getfield(h,'database()');
	database = database{1};
else
	try
		h = mym(['use ' database]);
	catch
		error(['Could not open database ''' database '''.']);
	end
end
