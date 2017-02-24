function DBUse(database)

%DBUse - Set current database.
%
% Set current database to store/retrieve processed data.
%
%  USAGE
%
%    DBUse(database)
%
%    database           database name
%
%  SEE
%
%    See also DBConnect.
%

% Copyright (C) 2007-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

try
	h = mym(['use ' database]);
catch
   error(['Could not open database ''' database '''.']);
end
