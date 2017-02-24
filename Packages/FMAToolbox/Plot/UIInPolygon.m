function in = UIInPolygon(X,Y)

%UIInPolygon - Find points in interactively defined polygon zone.
%
%  Mouse buttons:
%   left    =   add polygon point
%   right   =   remove last polygon point
%   middle  =   close polygon
%
%  USAGE
%
%    in = UIInPolygon(X,Y)
%
%    X,Y            polygon coordinates
%    in             logical vector indicating whether points at (X,Y) are inside
%
%  SEE
%
%    See also UISelect
%

% Copyright (C) 2009-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

[Xp,Yp,p] = UISelect;
in = inpolygon(X,Y,Xp,Yp);
hold on;
pts = plot(X(in),Y(in),'+r');
hold off;
pause(0.4);
delete(p);
delete(pts);
