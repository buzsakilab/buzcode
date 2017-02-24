function [v,comments,parameters,mfiles,code,date,user] = DBSelectVariables(query)

%DBSelectVariables - Select variables in current database.
%
% Select variables. Also get code used to generate them.
%
%  USAGE
%
%    x = DBSelectVariables(query)
%
%    query          optional list query (WHERE clause; see Example)
%
%  OUTPUT
%
%    x is a structure with the following fields:
%
%    eid            experiment ID (identifier string)
%    name           figure descriptive name (identifier string)
%    v              variable value
%    comments       comments
%    parameters     figure parameters
%    code           text of m-files used to generate the figure
%    date           date when the figure was saved
%    user           the user that saved the figure
%
%    Each field is a cell array.
%
%  EXAMPLE
%
%    Get all variables from experiment "experiment1", the name of which starts
%    with "raster" (for details, see an SQL manual):
%
%    v = DBSelectVariables('eid="experiment1" and name like "raster%"');

% Copyright (C) 2007-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

% Edit query
if nargin == 0,
	query = '';
end
query = strtrim(query);
query = regexprep(query,'^where','');
if ~isempty(query), query = [' where ' query]; end

% Query database
v = mym(['select v,eid,name,comments,parameters,mfiles,code,date,user from variables' query]);

if isempty(v),
	warning(['No variables match (' query ').']);
end

for i = 1:length(v.code),
	for j = 1:length(v.code{i}),
		code{i}{j} = char(v.code{i}{j})';
	end
end
v.code = code;
