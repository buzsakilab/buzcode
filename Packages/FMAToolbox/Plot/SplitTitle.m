function SplitTitle(h,text,width)

%SplitTitle - Split figure title over multiple lines
%
%  USAGE
%
%    SplitTitle(text,width)
%    SplitTitle(h,text,width)
%
%    h              optional handle (default = gca)
%    text           title text
%    width          optional line width (in chars, default = 10)
%

% Copyright (C) 2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help SplitTitle">SplitTitle</a>'' for details).');
end

% Default values
if nargin == 1,
	text = h;
	h = gca;
	width = 10;
elseif nargin == 2,
	if ishandle(h),
		width = 10;
	else
		width = text;
		text = h;
		h = gca;
	end
end

if ~isiscalar(width),
  error('Non-integer width (type ''help <a href="matlab:help SplitTitle">SplitTitle</a>'' for details).');
end

% Split text across N lines of fixed width
l = length(text);
nLines = ceil(l/width);
start = -width+1;
t = {};
for i = 1:nLines-1,
	start = (i-1)*width+1;
	t{i} = text(start:start+width-1);
end
t{end+1} = text(start+width:end);

% Set title
switch lower(get(h,'type')),
	case 'uipanel',
		set(h,'title',t);
	case 'axes',
		title(h,t);
end
