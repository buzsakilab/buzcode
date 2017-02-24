function csd = CSD(lfp)

%CSD - Compute current source density.
%
%  USAGE
%
%    csd = CSD(lfp)
%
%    lfp            local field potential samples
%
%  SEE
%
%    See also PlotCSD.

% Copyright (C) 2008-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CSD">CSD</a>'' for details).');
end

t = lfp(:,1);
y = lfp(:,2:end);
y = y - repmat(mean(y),length(t),1);
d = -diff(y,2,2);
csd = [t d];
