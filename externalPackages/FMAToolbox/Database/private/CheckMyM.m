%CheckMyM - Make sure mym is installed and functional.
%
% Connectivity to MySQL database server depends on mym. Make sure
% this is installed and functional.
%
%  USAGE
%
%    CheckMyM

% Copyright (C) 2007-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function CheckMyM

if isempty(which('mym')),
	error('FMAToolbox:CheckMyM:MissingMyM','This function requires the <a href="http://sourceforge.net/projects/mym/">mYm</a> toolbox by Y. Maret.\nTo compile it, you will need to install MySQL, zlib and the corresponding header files. Assuming\nthe libraries are in /usr/lib, and the headers are in /usr/include/mysql and /usr/include, type:\n\n  >> mex -v -I/usr/include/mysql -I/usr/include -L/usr/lib -lz -lmysqlclient mym.cpp');
end
