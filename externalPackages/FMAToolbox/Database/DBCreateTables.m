function DBCreateTables(database)

%DBCreateTables - Create tables in existing database.
%
% Create figure and variable tables in an existing database. If you
% have admin rights to the database server, use <a href="matlab:help DBCreate">DBCreate</a> instead.
% Otherwise ask the DBA to create an empty database with the appropriate
% access rights, and then use DBCreateTables to create the actual tables
% in the database.
%
%  USAGE
%
%    DBCreateTables(database)
%
%    database           database name
%
%  SEE
%
%    See also DBCreate.
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
  error('Incorrect number of parameters (type ''help <a href="matlab:help DBCreateTables">DBCreateTables</a>'' for details).');
end
if ~isstring_FMAT(database),
  error('Incorrect database name (type ''help <a href="matlab:help DBCreateTables">DBCreateTables</a>'' for details).');
end

h = mym(['use ' database]);

try
	h = mym('create table figures (eid varchar(50),name varchar(100),comments varchar(255),parameters varchar(50),mfiles mediumblob,code mediumblob,date varchar(50),user varchar(15),md5 varchar(32),gid smallint unsigned,fig mediumblob,png mediumblob,primary key (eid,name));');
	h = mym('create table variables (eid varchar(50),name varchar(100),comments varchar(255),parameters varchar(50),mfiles mediumblob,code mediumblob,date varchar(50),user varchar(15),md5 varchar(32),gid smallint unsigned,v longblob,primary key (eid,name));');
catch
   error(['Could not create tables in ''' database '''.']);
end
