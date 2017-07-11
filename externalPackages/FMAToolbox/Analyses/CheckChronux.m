%CheckChronux - Make sure chronux is installed and functional.
%
%  USAGE
%
%    CheckChronux

% Copyright (C) 2007-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function CheckChronux

if isempty(which('chronux')),
	error('This function requires the <a href="http://www.chronux.org">chronux</a> toolbox.');
end
