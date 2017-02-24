function x = DBGetVariables(query,varargin)

%DBGetVariables - Get all variables that match given criteria.
%
% Get variables and related information such as code used to generate them.
%
% This is mostly useful to get information about the variables (who stored
% them, when, which code was used, with what parameters, etc.) However, if
% you are only interested in the actual values for subsequent processing
% - which is the most frequent scenario -, it will be more appropriate to
% use <a href="DBGetValues">DBGetValues</a> instead, which provides a more useful output: the values
% will be arranged into a matrix or 3D array, rather than a cumbersome
% structure of cell arrays.
%
% There are exceptions, however. Some types of data, such as character strings
% of varying lengths or matrices of unequal sizes, cannot be grouped into a
% matrix or 3D array. In such cases, you will not be able to use <a href="DBGetValues">DBGetValues</a>,
% and will have to use DBGetVariables instead.
%
%  USAGE
%
%    x = DBGetVariables(query,<options>)
%
%    query          optional database query (WHERE clause; see Example)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'output'       'variables' to get the variables but not the information,
%                    i.e. eid, name, comments, parameters, etc. (default),
%                    'info' for the information but not the variables, 'full'
%                    for both, or 'keys' for eid and name only
%    =========================================================================
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
%    x = DBGetVariables('eid="experiment1" and name like "raster%"');
%
%  SEE
%
%    See also DBGetValues, DBAddVariable, DBGetFigures, DBDisplay.
%

% Copyright (C) 2007-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

% Default values
output = 'variables';

% Optional query provided?
if nargin == 0,
	query = '';
elseif isstring_FMAT(query,'output'),
	varargin = {query varargin{:}};
	query = '';
end

% Edit query
query = strtrim(query);
query = regexprep(query,'^where','');
if ~isempty(query), query = [' where ' query]; end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help DBGetVariables">DBGetVariables</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'output',
			output = lower(varargin{i+1});
			if ~isstring_FMAT(output,'variables','info','full','keys'),
				error('Incorrect value for property ''output'' (type ''help <a href="matlab:help DBGetVariables">DBGetVariables</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help DBGetVariables">DBGetVariables</a>'' for details).']);
	end
end

% Query database
switch output,
	case 'full',
		x = mym(['select v,eid,name,comments,parameters,mfiles,code,date,user from variables' query]);
	case 'variables',
		x = mym(['select v from variables' query]);
	case 'info',
		x = mym(['select eid,name,comments,parameters,mfiles,code,date,user from variables' query]);
	case 'keys',
		x = mym(['select eid,name from variables' query]);
end


% Make sure query results are not empty
if isempty(x),
	warning(['No variables match (' query ').']);
end

% Format code
if strcmp(output,'full') || strcmp(output,'info'),
	code = {};
	for i = 1:length(x.code),
		for j = 1:length(x.code{i}),
			code{i}{j} = char(x.code{i}{j})';
		end
	end
	x.code = code;
end