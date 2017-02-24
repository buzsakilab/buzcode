function [v,comments,parameters,mfiles,code,date,user] = DBSelectVariable(eid,name)

%DBSelectVariable - Select variable in current database.
%
% Select variable. Also get code used to generate it.
%
%  USAGE
%
%    x = DBSelectVariable(eid,name)
%
%    eid            experiment ID (identifier string)
%    name           variable descriptive name (identifier string)
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

% Copyright (C) 2007-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

% Check number of parameters
if nargin ~= 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help DBSelectVariable">DBSelectVariable</a>'' for details).');
end

v = mym(['select v,eid,name,comments,parameters,mfiles,code,date,user from variables where eid="' eid '" and name="' name '"']);
if isempty(v),
	warning(['Variable (' eid ',' name ') not found.']);
end

for i = 1:length(v.code{1}),
	code{i} = char(v.code{1}{i})';
end
v.code = code;
