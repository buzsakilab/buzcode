function databases = DBList(pattern)

%DBList - List existing databases.
%
%  USAGE
%
%    databases = DBList(pattern)
%
%    pattern        optional name search pattern (see below)
%
%  EXAMPLES
%
%    % List all databases
%    databases = DBList;
%
%    % List all databases starting with 'LFP'
%    databases = DBList('LFP.*');
%
%  SEE
%
%    See also DBConnect.
%

% Copyright (C) 2012-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

try
	h = mym('show databases;');
catch
	error(['Could not list databases (check DB server connection).']);
end
databases = h.Database;

if nargin >= 1,
	databases = regexp(databases,pattern,'match');
	databases = vertcat(databases{:});
end
