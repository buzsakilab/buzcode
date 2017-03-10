function DBAddFigure(f,eid,name,comments,parameters,mfiles,varargin)

%DBAddFigure - Save figure in current database.
%
% Save figure as MATLAB fig and PNG image. Also save code used to generate it
% (although not required, installing <a href="http://sites.google.com/site/oliverwoodford/software/export_fig">export_fig</a> will yield better PNG images).
%
%  USAGE
%
%    DBAddFigure(f,eid,name,comments,parameters,mfiles,<options>)
%
%    f              figure handle
%    eid            experiment ID (identifier string)
%    name           descriptive name
%    comments       comment string
%    parameters     figure parameters
%    mfiles         cell array of m-file names (used to generate the figure)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'onlyPNG'     skip fig and only store png capture (use this to avoid
%                   storing very large figures into the database,
%                   default = 'off')
%    =========================================================================
%
%  NOTE
%
%    The pair (eid,name) must uniquely define the figure (they are used as
%    primary keys in the database).
%
%  SEE
%
%    See also DBAddVariable, DBGetFigures, DBExportGallery.
%

% Copyright (C) 2007-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global dbUser;

% Make sure MyM is installed and functional
CheckMyM;

% Check number of parameters
if nargin < 6,
	error('Incorrect number of parameters (type ''help <a href="matlab:help DBAddFigure">DBAddFigure</a>'' for details).');
end

% Default values
onlyPNG = 'off';

% Check parameters
if ~ishandle(f),
	error('Incorrect figure handle (type ''help <a href="matlab:help DBAddFigure">DBAddFigure</a>'' for details).');
end
if ~ischar(eid),
	error('Incorrect EID string (type ''help <a href="matlab:help DBAddFigure">DBAddFigure</a>'' for details).');
end
if ~ischar(name),
	error('Incorrect figure name (type ''help <a href="matlab:help DBAddFigure">DBAddFigure</a>'' for details).');
end
if ~isempty(strfind(eid,'/')) || ~isempty(strfind(name,'/')),
	error('EIDs and names should not contain ''/'' (type ''help <a href="matlab:help DBAddFigure">DBAddFigure</a>'' for details).');
end
if ~ischar(comments),
	error('Incorrect comment string (type ''help <a href="matlab:help DBAddFigure">DBAddFigure</a>'' for details).');
end
if ~ischar(parameters),
	error('Incorrect figure parameters (type ''help <a href="matlab:help DBAddFigure">DBAddFigure</a>'' for details).');
end
if ~ischar(mfiles) && ~iscell(mfiles),
	error('Incorrect list of m-files (type ''help <a href="matlab:help DBAddFigure">DBAddFigure</a>'' for details).');
end
if ischar(mfiles),
	mfiles = {mfiles};
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+firstIndex) ' is not a property (type ''help <a href="matlab:help DBAddFigure">DBAddFigure</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'onlypng',
			onlyPNG = lower(varargin{i+1});
			if ~isstring_FMAT(onlyPNG,'on','off'),
				error('Incorrect value for property ''onlypng'' (type ''help <a href="matlab:help DBAddFigure">DBAddFigure</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help DBAddFigure">DBAddFigure</a>'' for details).']);
	end
end

onlyPNG = strcmp(onlyPNG,'on');

basename = tempname;
figName = [basename '.fig'];
pngName = [basename '.png'];
if ~onlyPNG,
	try
		% Create temporary fig file
		saveas(f,figName);
	catch
		error('Could not create temporary file for fig figure.');
	end
end
try
	% Create temporary png file
	if ~isempty(which('export_fig')),
		export_fig(f,pngName,'-png','-a1');
	else
		saveas(f,pngName);
	end
catch
	error('Could not create temporary file for png figure.');
end

% M-Files
for i = 1:length(mfiles),
	mfileName = which(mfiles{i});
	if isempty(mfileName),
		error(['M-File for function ''' mfiles{i} ''' not found.']);
	end
	fid = fopen(mfileName);
	code{i} = fread(fid);
	fclose(fid);
end

% Date and md5
d = datestr(now);
md5 = CalcMD5(pngName);

% Execute SQL query command
if onlyPNG,
	h = mym(['insert into figures (eid,name,comments,parameters,mfiles,code,date,user,md5,png) values ("{S}","{S}","{S}","{S}","{M}","{M}","{S}","{S}","{Si}","{uF}")'],eid,name,comments,parameters,mfiles,code,d,dbUser,md5,pngName);
else
	h = mym(['insert into figures (eid,name,comments,parameters,mfiles,code,date,user,md5,fig,png) values ("{S}","{S}","{S}","{S}","{M}","{M}","{S}","{S}","{Si}","{F}","{uF}")'],eid,name,comments,parameters,mfiles,code,d,dbUser,md5,figName,pngName);
end

% Delete temporary files
if ~onlyPNG,
	delete(figName);
end
delete(pngName);
