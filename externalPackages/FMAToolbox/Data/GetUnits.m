function units = GetUnits(group)

%GetUnits - Get list of units.
%
%  USAGE
%
%    units = GetUnits(group)
%
%    group          optional electrode group (default = all groups)
%
%  EXAMPLES
%
%    all = GetUnits;        % list all units
%    units = GetUnits(2);   % list all units on electrode group 2


% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA;
if isempty(DATA),
	error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end

% Exclude clusters 0 and 1 (artefacts and noise)
units = DATA.spikes(:,[2 3]);
units = unique(units,'rows');
units = units(units(:,2)~=0&units(:,2)~=1,:);

if nargin > 0,
	units = units(units(:,1)==group,:);
end
