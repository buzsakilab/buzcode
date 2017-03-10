function database = DBCreate(database)

%DBCreate - Create database.
%
% Create database with figure and variable tables. Use this function
% if you have admin rights to the database server. Otherwise ask
% the DBA to create an empty database with the appropriate access
% rights, and then use <a href="matlab:help DBCreateTables">DBCreateTables</a> to create the actual tables in
% the database.
%
%  USAGE
%
%    name = DBCreate(database)
%
%    database       database name (see NOTE below)
%
%  NOTE
%
%    Database names can include wildcards to indicate current date and time:
%
%      %y    year
%      %m    month
%      %d    day
%      %t    time
%
%  EXAMPLE
%
%    DBCreate('TestData_%y%m%d');
%
%  SEE
%
%    See also DBCreateTables, DBRemove.
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
  error('Incorrect number of parameters (type ''help <a href="matlab:help DBCreate">DBCreate</a>'' for details).');
end
if ~isstring_FMAT(database),
  error('Incorrect database name (type ''help <a href="matlab:help DBCreate">DBCreate</a>'' for details).');
end

% Insert date
database = InsertDate(database);

% Create database
try
	h = mym(['create database ' database]);
catch
   error('FMAToolbox:DBCreate:cannotCreateDB',['Could not create database ''' database '''.\nThis database might already exist or its name may contain incorrect characters (e.g. ''-'').\nAlternatively, you may not have the appropriate access rights, in which case you may want\nto ask your DBA to create an empty database for you, and then use <a href="matlab:help DBCreateTables">DBCreateTables</a> to create\nthe tables in the database.']);
end

h = mym(['use ' database]);

% Create tables
try
	h = mym('create table figures (eid varchar(50),name varchar(100),comments varchar(255),parameters varchar(50),mfiles mediumblob,code mediumblob,date varchar(50),user varchar(15),md5 varchar(32),gid smallint unsigned,fig mediumblob,png mediumblob,primary key (eid,name));');
	h = mym('create table variables (eid varchar(50),name varchar(100),comments varchar(255),parameters varchar(50),mfiles mediumblob,code mediumblob,date varchar(50),user varchar(15),md5 varchar(32),gid smallint unsigned,v longblob,primary key (eid,name));');
catch
   error(['Could not create tables in ''' database '''.']);
end
