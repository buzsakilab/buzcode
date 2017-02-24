function f = DBSelectFigures(query)

%DBSelectFigures - Select figures in current database.
%
% Select figures and open as MATLAB fig. Also get code used to generate them.
%
%  USAGE
%
%    f = DBSelectFigures(query)
%
%    query          optional figure list query (WHERE clause; see Example)
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
%      DBSelectFigures('eid="experiment1" and name like "raster%"');
%
%  SEE
%
%    See also DBSelectFigure.
%

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
f = mym(['select fig,eid,name,comments,parameters,mfiles,code,date,user from figures' query]);
if isempty(f),
	warning(['No figures match (' query ').']);
end

for i = 1:length(f.code),
	for j = 1:length(f.code{i}),
		code{i}{j} = char(f.code{i}{j})';
	end
end
f.code = code;

% Create temporary fig files, open them and delete them
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