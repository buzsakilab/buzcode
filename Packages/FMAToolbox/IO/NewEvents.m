function events = NewEvents(times,description)

%NewEvents - Create events structure.
%
%  USAGE
%
%    events = NewEvents(times,description)
%
%    times               event timestamps
%    description         (common) event description

% Copyright (C) 2008-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

events.time = times;
events.description = cell(size(times));
[events.description{:}] = deal(description);
