function f = DBSelectFigure(eid,name)

%DBSelectFigure - Select figure in current database.
%
% Select figure and open as MATLAB fig. Also get code used to generate it.
%
%  USAGE
%
%    f = DBSelectFigure(eid,name)
%
%    eid            experiment ID (identifier string)
%    name           figure descriptive name (identifier string)
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
%  SEE
%
%    See also DBSelectFigures.
%

% Copyright (C) 2007-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

% Check number of parameters
if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help DBSelectFigure">DBSelectFigure</a>'' for details).');
end

f = mym(['select fig,comments,parameters,mfiles,code,date,user from figures where eid="' eid '" and name="' name '"']);
if isempty(f.fig),
	error(['Figure (' eid ',' name ') not found.']);
end
if isempty(f.fig{1}),
	error(['Figure (' eid ',' name ') was stored as PNG only.']);
end

for i = 1:length(f.code{1}),
	code{i} = char(f.code{1}{i})';
end
f.code = code;

f.eid = {eid};
f.name = {name};

% Create temporary fig file, open it and delete it
basename = tempname;
figName = [basename '.fig'];
file = fopen(figName,'wb');
if file == -1,
	error('Could not create temporary file for figure.');
end
fwrite(file,f.fig{1});
fclose(file);
f.fig = {openfig(figName)};
delete(figName);
set(f.fig{1},'visible','on');
