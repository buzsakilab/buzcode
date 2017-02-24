function out = HTML(in,format)

%HTML - Create HTML formatted string.
%
%  Important note: Matlab UITables are buggy! Therefore, HTML formatted table
%  cells will not be rendered satisfactorily (e.g. you will experience problems
%  with alignment and cell width). This function does not try to compensate for
%  these bugs.
%
%  USAGE
%
%    out = HTML(in,format)
%
%    in             string, number or cell array to format
%    format         color, weight and angle information
%
%  NOTE
%
%    The format string is a concatenation of one or more of the following items:
%
%      B           boldface
%      I           italics
%      #nnnnnn     HTML color code (e.g. #f1c27a)
%
%  EXAMPLES
%
%    % Create a table figure with two items,
%    % the first one in bright red bold characters
%
%    m = { HTML('Total','B#ff0000'),200 };
%    TableFigure('Results',m);
%
%    % Create a table figure with two items,
%    % both in bright red bold characters
%
%    m = HTML({'Total',200},'B#ff0000');
%    TableFigure('Results',m);
%
%  SEE
%
%    See also TableFigure.

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check parameters
if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help HTML">HTML</a>'' for details).');
end
if ~iscell(in),
	in = {in};
	inputIsCell = 0;
else
	inputIsCell = 1;
end
for i = 1:length(in(:)),
	if ~ischar(in{i}) && ~isscalar(in{i}),
		error('Incorrect string or number (type ''help <a href="matlab:help HTML">HTML</a>'' for details).');
	end
end
if ~ischar(format),
  error('Incorrect format (type ''help <a href="matlab:help HTML">HTML</a>'' for details).');
end

format = lower(format);
out = in;

% Parse format string
r = regexp(format,'([bi]*)((#[0-9a-f]{6})*)([bi]*)','tokens');
if isempty(r), return; end
r = r{1};
if ~strcmp(format,[r{:}]),
  error('Incorrect format (type ''help <a href="matlab:help HTML">HTML</a>'' for details).');
end

global color style;

% Color
c = r{2};
if isempty(c),
	color = '';
else
	color = ['color: ' c ';'];
end

% Style
s = [r{[1 3]}];
style = '';
if ~isempty(s),
	if any(s=='i'),
		style = 'font-angle: italic;';
	end
	if any(s=='b'),
		style = [style 'font-weight: bold;'];
	end
end

if isempty([color style]), return; end

% Output HTML string
out = cellfun(@FormatCell,in,'uniformoutput',0);

if ~inputIsCell,
	out = out{:};
end


function y = FormatCell(x)

global color style;

if isnumeric(x),
	x = num2str(x);
	align = 'text-align:right;';
else
	align = '';
end
y = ['<html><span style="' color style align '">' x '</span></html>'];
