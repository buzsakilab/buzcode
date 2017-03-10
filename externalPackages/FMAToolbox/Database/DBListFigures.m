function DBListFigures(query)

%DBListFigures - List figures matching certain criteria.
%
%  USAGE
%
%    DBListFigures(query)
%
%    query          optional database query (WHERE clause)
%
%  EXAMPLE
%
%    DBListFigures('name="Spectrogram"');
%
%  SEE
%
%    See also DBGetVariables, DBGetFigures.
%

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

if nargin < 1,
	query = '';
end

DBDisplay(DBGetFigures(query,'output','keys'));
