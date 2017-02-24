function [X,Y,polygon] = UISelect

%UISelect - Interactively select polygon zone in an existing plot.
%
%  USAGE
%
%    [X,Y,p] = UISelect
%
%    X,Y            polygon coordinates
%    p              handle for polygon object in figure
%
%  SEE
%
%    See also UIInPolygon
%

% Copyright (C) 2009-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

hold on;
X = [];
Y = [];
polygon = [];
last = false;

while ~last,
	% Get next mouse click
	[x,y,button] = ginput(1);
	switch button,
		case 1,
			% Left button adds this point
			X(end+1,1) = x;
			Y(end+1,1) = y;
		case 2,
			% Middle button closes the selection
			last = true;
			if ~isempty(X),
				X(end+1,1) = X(1);
				Y(end+1,1) = Y(1);
			end
		case 3,
			% Right button cancels last point
			if ~isempty(X),
				X(end) = [];
				Y(end) = [];
			end
	end
	% Update existing polygon
	if isempty(polygon),
		polygon = plot(X,Y);
	else
		warning('off','MATLAB:hg:line:XDataAndYDataLengthsMustBeEqual');
		set(polygon,'XData',X);
		set(polygon,'YData',Y);
		warning('on','MATLAB:hg:line:XDataAndYDataLengthsMustBeEqual');
	end
end

hold off;
