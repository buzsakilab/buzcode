function DBExportGallery(query,name,varargin)

%DBExportGallery - Create figure gallery from current database.
%
% Select figures and create an HTML gallery (although not required, installing
% <a href="http://sites.google.com/site/oliverwoodford/software/export_fig">export_fig</a> will yield better PNG images).
%
%  USAGE
%
%    DBExportGallery(query,name,<options>)
%
%    query          figure list query (WHERE clause; see Example)
%    name           gallery name
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'size'        horizontal and vertical thumbnail size
%                   (default = [240 240])
%     'nColumns'    number of columns (default = 3)
%     'info'        display detailed info with thumbnails (default = 'off')
%     'code'        create generation code pages (default = 'off')
%    =========================================================================
%
%  EXAMPLES
%
%    % Export entire gallery with extensive information
%    DBExportGallery('','Complete Gallery','info','on');
%
%    % Export only figures for the session named 'SESSION1'
%    DBExportGallery('eid="SESSION1"','Partial Gallery');
%
%    % Export figures for all sessions with names starting with 'SLEEP'
%    DBExportGallery('eid like "SLEEP%"','Sleep Gallery');

% Copyright (C) 2007-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

% Default values
hSize = 240;
vSize = 240;
nColumns = 3;
showInfo = 'off';
showCode = 'off';

% Check number of parameters
if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help DBExportGallery">DBExportGallery</a>'' for details).');
end

% Parse parameter list
for j = 1:2:length(varargin),
	if ~ischar(varargin{j}),
		error(['Parameter ' num2str(j+7) ' is not a property (type ''help <a href="matlab:help DBExportGallery">DBExportGallery</a>'' for details).']);
	end
	switch(lower(varargin{j})),
		case 'size',
			sizes = varargin{j+1};
			if ~isivector(sizes,'>0','#2'),
				error('Incorrect value for ''size'' (type ''help <a href="matlab:help DBExportGallery">DBExportGallery</a>'' for details).');
			end
			hSize = sizes(1);
			vSize = sizes(2);

		case 'ncolumns',
			nColumns = varargin{j+1};
			if ~isiscalar(nColumns,'>0'),
				error('Incorrect value for ''ncolumns'' (type ''help <a href="matlab:help DBExportGallery">DBExportGallery</a>'' for details).');
			end

		case 'info',
			showInfo = lower(varargin{j+1});
			if ~isstring_FMAT(showInfo,'on','off'),
				error('Incorrect value for ''info'' (type ''help <a href="matlab:help DBExportGallery">DBExportGallery</a>'' for details).');
			end

		case 'code',
			showCode = lower(varargin{j+1});
			if ~isstring_FMAT(showCode,'on','off'),
				error('Incorrect value for ''code'' (type ''help <a href="matlab:help DBExportGallery">DBExportGallery</a>'' for details).');
			end

		otherwise,
			error(['Unknown property ''' num2str(varargin{j}) ''' (type ''help <a href="matlab:help DBExportGallery">DBExportGallery</a>'' for details).']);

	end
end

showInfo = strcmp(showInfo,'on');
showCode = strcmp(showCode,'on');

% Edit query
query = strtrim(query);
query = regexprep(query,'^where','');
if ~isempty(query), query = [' where ' query]; end

% Query database
f = mym(['select png,eid,name,parameters,comments,mfiles,code,date,user from figures' query]);
if isempty(f),
	warning(['No figures match (' query ').']);
end

for i = 1:length(f.code),
	for j = 1:length(f.code{i}),
		code{i}{j} = char(f.code{i}{j})';
	end
end
f.code = code;

% Create root directory
try
	mkdir(name);
	mkdir([name '/images']);
	mkdir([name '/thumbs']);
	if showCode,
		mkdir([name '/code']);
	end
catch
	error('Cannot create gallery (check file access permissions).');
end
[path,name] = fileparts(name);

nFigures = length(f.png);

% Export png images + thumbnails (and code pages if required)
for i = 1:nFigures,
	% png image
	figureName{i} = [f.eid{i} '-' f.name{i} '.png'];
	file = fopen([path name '/images/' figureName{i}],'wb');
	if file == -1,
		error(['Could not create figure ' figureName{i} '.']);
	end
	fwrite(file,f.png{i});
	fclose(file);

	% Thumbnail
	pngData = imread([path name '/images/' figureName{i}],'png');
   [height,width,nz] = size(pngData);
   hScale = width/hSize;
   vScale = height/vSize;
	x = floor(1:hScale:width);
	y = floor(1:vScale:height);
   thumbnail = pngData(y,x,:);
   imwrite(thumbnail,[path name '/thumbs/' figureName{i}],'png');

	% Code
	if showCode && ~isempty(f.code{i}),
		file = fopen([path name '/code/' figureName{i}(1:end-4) '.m'],'w');
		if file == -1,
			error('Could not create code output.');
		end
		nFunctions = size(f.code{i},2);
		for j = 1:nFunctions,
			if j > 1, fprintf(file,'\n\n'); end
			fprintf(file,'%% =======================================================================================================\n');
			fprintf(file,['%%     ' f.mfiles{i}{j} '\n']);
			fprintf(file,'%% =======================================================================================================\n\n\n');
			fprintf(file,'%s',f.code{i}{j});
		end
%  		if isempty(j), j = 1; end
%  		fprintf(file,'%s',f.code{i}{j});
		fclose(file);
	end
end

% Output HTML code
% Header...
fid = fopen([path name '/index.html'],'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8" ?>\n');
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">\n');
fprintf(fid,'<html xmlns="http://www.w3.org/1999/xhtml">\n');
fprintf(fid,' <head>\n');
fprintf(fid,'  <title>%s</title>\n',name);
fprintf(fid,'  <meta http-equiv="content-type" content="text/html; charset=UTF-8"/>\n');
fprintf(fid,'  <meta name="GENERATOR" content="FMAToolbox"/>\n');
fprintf(fid,'  <style type="text/css">\n');
fprintf(fid,'   BODY                {color: #dddddd; background: #333333;margin: 1%%;\n');
fprintf(fid,'                        font-family: sans-serif; font-size: 10pt; }\n');
fprintf(fid,'   A:LINK              {color: #dddddd;}\n');
fprintf(fid,'   A:VISITED           {color: #cccccc;}\n');
fprintf(fid,'   H1                  {color: #eeeeee;}\n');
fprintf(fid,'   TABLE.enclosing     {margin-left: auto; margin-right: auto}\n');
fprintf(fid,'   TD.enclosing        {text-align: center; vertical-align:top; color: #dddddd; padding: 10px;background: #555555;border: 5px solid #333333;width: %dpx}\n',hSize*1.5);
fprintf(fid,'   TABLE.enclosed      {color: #dddddd; padding: 10px 0px 0px 0px;background: #555555;width: %dpx}\n',hSize*1.5-20);
fprintf(fid,'   TD.left             {text-align: right; color: #dddddd; padding: 0.1em 1em;background: #606060;border: 0px}\n');
fprintf(fid,'   TD.right            {text-align: left; color: #dddddd; padding: 0.1em 1em;background: #777777;border: 0px}\n');
fprintf(fid,'   IMG                 {border: 1px solid #dddddd;}\n');
fprintf(fid,'  </style>\n');
fprintf(fid,' </head>\n');
fprintf(fid,' <body>\n');
fprintf(fid,'  <h1>%s</h1>\n',name);
fprintf(fid,'  <p>%d images (%s)<br/></p>\n',nFigures,datestr(now));
fprintf(fid,'  <hr/>\n');
fprintf(fid,'  <table class="enclosing">\n');

% ... figure table...
nLines = ceil(nFigures/nColumns);
i = 1;
for line = 1:nLines,
	fprintf(fid,'   <tr>\n');
	for column = 1:nColumns,
		if i > nFigures, break; end
		fprintf(fid,'    <td class="enclosing">\n');
		fprintf(fid,'     <a href="images/%s"><img src="thumbs/%s" width="%d" height="%d"></a>\n',figureName{i},figureName{i},hSize,vSize);
		fprintf(fid,'     <table class="enclosed">\n');
		fprintf(fid,'      <tr>\n');
		fprintf(fid,'       <td class="left">\n');
		fprintf(fid,'        EID\n');
		fprintf(fid,'       </td>\n');
		fprintf(fid,'       <td class="right">\n');
		fprintf(fid,'        <div>%s</div>\n',f.eid{i});
		fprintf(fid,'       </td>\n');
		fprintf(fid,'      </tr>\n');
		fprintf(fid,'      <tr>\n');
		fprintf(fid,'       <td class="left">\n');
		fprintf(fid,'        Name\n');
		fprintf(fid,'       </td>\n');
		fprintf(fid,'       <td class="right">\n');
		fprintf(fid,'        <div>%s</div>\n',f.name{i});
		fprintf(fid,'       </td>\n');
		fprintf(fid,'      </tr>\n');
		if showInfo,
			fprintf(fid,'      <tr>\n');
			fprintf(fid,'       <td class="left">\n');
			fprintf(fid,'        Comments\n');
			fprintf(fid,'       </td>\n');
			fprintf(fid,'       <td class="right">\n');
			fprintf(fid,'        <div>%s</div>\n',f.comments{i});
			fprintf(fid,'       </td>\n');
			fprintf(fid,'      </tr>\n');
			fprintf(fid,'      <tr>\n');
			fprintf(fid,'       <td class="left">\n');
			fprintf(fid,'        Parameters\n');
			fprintf(fid,'       </td>\n');
			fprintf(fid,'       <td class="right">\n');
			fprintf(fid,'        <div>%s</div>\n',f.parameters{i});
			fprintf(fid,'       </td>\n');
			fprintf(fid,'      </tr>\n');
		end
		if showCode && ~isempty(f.code{i}),
			codeName = [figureName{i}(1:end-4) '.m'];
			fprintf(fid,'      <tr>\n');
			fprintf(fid,'       <td class="left">\n');
			fprintf(fid,'        Code\n');
			fprintf(fid,'       </td>\n');
			fprintf(fid,'       <td class="right">\n');
			fprintf(fid,'        <a href="code/%s">[M-File]</a>\n',codeName);
			fprintf(fid,'       </td>\n');
			fprintf(fid,'      </tr>\n');
		end
		if showInfo,
			fprintf(fid,'      <tr>\n');
			fprintf(fid,'       <td class="left">\n');
			fprintf(fid,'        Date\n');
			fprintf(fid,'       </td>\n');
			fprintf(fid,'       <td class="right">\n');
			fprintf(fid,'        <div>%s</div>\n',f.date{i});
			fprintf(fid,'       </td>\n');
			fprintf(fid,'      </tr>\n');
			fprintf(fid,'      <tr>\n');
			fprintf(fid,'       <td class="left">\n');
			fprintf(fid,'        User\n');
			fprintf(fid,'       </td>\n');
			fprintf(fid,'       <td class="right">\n');
			fprintf(fid,'        <div>%s</div>\n',f.user{i});
			fprintf(fid,'       </td>\n');
			fprintf(fid,'      </tr>\n');
		end
		fprintf(fid,'     </table>\n');
		fprintf(fid,'    </td>\n');
		i = i+1;
	end
	fprintf(fid,'   </tr>\n');
end

% ... and footer
fprintf(fid,'  </table>\n');
fprintf(fid,' </body>\n');
fprintf(fid,'</html>\n');
