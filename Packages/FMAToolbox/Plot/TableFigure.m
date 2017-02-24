function [f,t] = TableFigure(title,data,columns)

%TableFigure - Create table figure.
%
%  USAGE
%
%    [f,t] = TableFigure(title,data,columns)
%
%    title          table title
%    data           cell array of strings to populate the table
%    column         optional cell vector of column names
%
%  OUTPUT
%
%    f              figure handle
%    t              table object
%
%  SEE
%
%    See also HTML for table cell format.

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help TableFigure">TableFigure</a>'' for details).');
end
if nargin < 3,
	columns = cell(1,size(data,2));
end

if size(columns,2) ~= size(data,2),
	error('Incompatible data and column sizes (type ''help <a href="matlab:help TableFigure">TableFigure</a>'' for details).');
end

titleHeight = 20;

f = figure;
pos = get(f,'position');
uicontrol('style','text','position',[0 pos(4)-titleHeight pos(3) titleHeight],'string',title);
p = uipanel('position',[0 0 pos(3) pos(4)-titleHeight]);
t = uitable('ColumnName',columns,'position',[0 0 pos(3) pos(4)-titleHeight],'data',data,'parent',p);

% Adjust column widths to contents

% Determine maximum width in each column
widths = max(cellfun(@(x) length(num2str(x)),[data;columns]));
% Convert to pixels
a1 = get(t,'extent');
a1 = a1(3);
set(t,'units','points');
a2 = get(t,'extent');
a2 = a2(3);
set(t,'units','pixels');
widths = round(1.2*widths/a1*a2*get(t,'fontsize'));
widths = num2cell(widths);
% Which columns contain numbers?
numbers = logical(sum(cellfun(@isnumeric,data)));
% Leave these columns in 'auto' mode
widths(numbers) = {'auto'};
% Adjust
set(t,'ColumnWidth',widths);