function V = LinearVelocity(X,smooth)

%LinearVelocity - Compute instantaneous linear velocity.
%
% Compute linear velocity for a time-varying vector X.
%
%  USAGE
%
%    V = LinearVelocity(X,smooth)
%
%    X              the time-varying vector given as a list of triplets [t x y]
%    smooth         optional standard deviation for Gaussian kernel used for
%                   differentiating, measured in number of samples
%                   (default = no smoothing)
%
%  SEE
%
%    See also AngularVelocity.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help LinearVelocity">LinearVelocity</a>'' for details).');
end
if nargin >= 2,
	if ~isdscalar(smooth,'>=0'),
		error('Incorrect smoothing stdev (type ''help <a href="matlab:help LinearVelocity">LinearVelocity</a>'' for details).');
	end
else
	smooth = 0;
end

try
DX = Diff(X,'smooth',smooth);
catch
DX(:,1) = Diff(X(:,1),'smooth',smooth);
DX(:,2) = Diff(X(:,2),'smooth',smooth);
DX(:,3) = Diff(X(:,3),'smooth',smooth);
end

Y = DX(:,2:3).*DX(:,2:3);
N = sqrt(Y(:,1)+Y(:,2));
V = [X(:,1) N];
