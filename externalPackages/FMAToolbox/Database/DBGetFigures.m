function f = DBGetFigures(query,varargin)

%DBGetFigures - Get all figures that match given criteria.
%
% Open figures and get related information such as code used to generate them.
%
%  USAGE
%
%    f = DBGetFigures(query,<options>)
%
%    query          optional figure list query (WHERE clause; see Example)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'output'       'figures' to get the figures but not the information,
%                    i.e. eid, name, comments, parameters, etc. (default),
%                    'info' for the information but not the figures, 'full'
%                    for both, or 'keys' for eid and name only
%    =========================================================================
%
%  OUTPUT
%
%    f is a structure with the following fields:
%
%    eid            experiment ID (identifier string)
%    name           figure descriptive name (identifier string)
%    fig            figure handle
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
%    Get all figures from experiment "experiment1", the name of which starts
%    with "raster" (for details, see an SQL manual):
%
%    [f,comments,parameters,code,date,user] = ...
%      DBGetFigures('eid="experiment1" and name like "raster%"');
%
%  SEE
%
%    See also DBAddFigure, DBGetVariables, DBDisplay.
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
output = 'figures';

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
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help DBGetFigures">DBGetFigures</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'output',
			output = lower(varargin{i+1});
			if ~isstring_FMAT(output,'figures','info','full','keys'),
				error('Incorrect value for property ''output'' (type ''help <a href="matlab:help DBGetFigures">DBGetFigures</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help DBGetFigures">DBGetFigures</a>'' for details).']);
	end
end

% Query database
switch output,
	case 'full',
		f = mym(['select fig,eid,name,comments,parameters,mfiles,code,date,user from figures' query]);
	case 'figures',
		f = mym(['select fig from figures' query]);
	case 'info',
		f = mym(['select eid,name,comments,parameters,mfiles,code,date,user from figures' query]);
	case 'keys',
		f = mym(['select eid,name from figures' query]);
end

% Make sure query results are not empty
if isempty(f),
	warning(['No figures match (' query ').']);
end

% Create temporary fig files, open them and delete them
if strcmp(output,'full') || strcmp(output,'figures'),
	for i = 1:length(f.fig),
		if isempty(f.fig{i}),
			warning(['Figure (' f.eid{i} ',' f.name{i} ') was stored as PNG only.']);
			f.fig{i} = [];
		else
			basename = tempname;
			figName = [basename '.fig'];
			file = fopen(figName,'wb');
			if file == -1,
				error(['Could not create temporary file for figure (' f.eid{i} ',' f.name{i} ').']);
			end
			fwrite(file,f.fig{i});
			fclose(file);
			f.fig{i} = openfig(figName);
			delete(figName);
			set(f.fig{i},'visible','on');
		end
	end
end

% Format code
if strcmp(output,'full') || strcmp(output,'info'),
	code = {};
	for i = 1:length(f.code),
		for j = 1:length(f.code{i}),
			code{i}{j} = char(f.code{i}{j})';
		end
	end
	f.code = code;
end
