%CheckChronux - Make sure chronux is installed and functional.
%
%  USAGE
%
%    CheckChronux
%       -or-
%    CheckChronux(functionname) where function name is a string of the
%                               chronux function you want to check

% Copyright (C) 2007-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function CheckChronux(functionname)

if isempty(which('chronux')),
	error('This function requires the <a href="http://www.chronux.org">chronux</a> toolbox.');
end

if exist('functionname','var')
    if isempty(which(functionname)),
        error(['Chronux function: ',functionname,' does not exist :(']);
    end
end

end
