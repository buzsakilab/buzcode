function list = DBListVariables

%DBListVariables - List variables in current database.
%
%  USAGE
%
%    list = DBListVariables
%
%  SEE
%
%    See also DBListFigures.
%

% Copyright (C) 2007-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

list = mym('select eid,name,comments,parameters,date,user from variables');
