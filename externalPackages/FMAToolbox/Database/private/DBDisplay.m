function DBDisplay(results)

%DBDisplay - Display results of a database query.
%
%  USAGE
%
%    DBDisplay(results)
%
%    results        output of <a href="matlab:help DBGetFigures">DBGetFigures</a> or <a href="matlab:help DBGetFigures">DBGetVariables</a>
%
%  EXAMPLE
%
%    DBDisplay(DBGetVariables('name="Spectrogram",'output','info'));
%
%  SEE
%
%    See also DBGetVariables, DBGetFigures.
%

% Copyright (C) 2007-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help DBDisplay">DBDisplay</a>'' for details).');
end

if ~isfield(results,'eid') || ~isfield(results,'name'),
	error('Incorrect query results (type ''help <a href="matlab:help DBDisplay">DBDisplay</a>'' for details).');
end

% Display

% Column widths
title1 = 'EID';
title2 = 'NAME';
if isempty(results.eid),
	width1 = 0;
	width2 = 0;
else
	width1 = max(cellfun(@length,results.eid));
	width2 = max(cellfun(@length,results.name));
end
width1 = max([width1 length(title1)]);
width2 = max([width2 length(title2)]);
spacing = '     ';
total = width1+width2+3*length(spacing);
% Print header
disp(' ');
disp(repmat('=',1,total));
disp([spacing sprintf(['%-' int2str(width1) 's%s%-' int2str(width2) 's'],'EID',spacing,'NAME')]);
disp(repmat('-',1,total'));
% Print rows
for i = 1:length(results.eid),
	disp([spacing sprintf(['%-' int2str(width1) 's%s%-' int2str(width2) 's'],results.eid{i},spacing,results.name{i})]);
end
% Print footer
disp(repmat('=',1,total'));
disp(' ');
