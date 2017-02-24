function out = InsertDate(in)

%InsertDate - Insert date into string.
%
% Replaces %y, %m, %d and %t with current year, month, day and time.
%
%  USAGE
%
%    out = InsertDate(in)
%
%    in             string to modify

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

n = clock;
year = int2str(n(1));
month = int2zstr(n(2),2);
day = int2zstr(n(3),2);
time = [int2zstr(n(4)) '' int2zstr(n(5))];

out = in;
out = strrep(out,'%y',year);
out = strrep(out,'%m',month);
out = strrep(out,'%d',day);
out = strrep(out,'%t',time);
