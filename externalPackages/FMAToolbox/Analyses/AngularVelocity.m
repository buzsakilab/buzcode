function W = AngularVelocity(X,smooth)

%AngularVelocity - Compute instantaneous angular velocity.
%
% Compute angular velocity for a time-varying vector X. This is computed as
% the vector product of X and its derivative DX: W(t) = X(t) ^ DX(t)
%
%  USAGE
%
%    omega = AngularVelocity(X,smooth)
%
%    X              the time-varying vector given as a list of triplets [t x y]
%    smooth         optional standard deviation for Gaussian kernel used for
%                   differentiating, measured in number of samples
%                   (default = no smoothing)
%
%  NOTE
%
%    Angular velocities are returned in radians/s.
%
%  SEE
%
%    See also LinearVelocity.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help AngularVelocity">AngularVelocity</a>'' for details).');
end
if nargin >= 2,
	if ~isdscalar(smooth,'>=0'),
		error('Incorrect smoothing stdev (type ''help <a href="matlab:help AngularVelocity">AngularVelocity</a>'' for details).');
	end
else
	smooth = 0;
end

% Norm X
Y = X(:,2:3).*X(:,2:3);
N = sqrt(Y(:,1)+Y(:,2));
X = [X(:,1) X(:,2)./N X(:,3)./N];

% Differentiate
DX = Diff(X,'smooth',smooth);

% Compute vectorial product

W = [X(:,1) X(:,2).*DX(:,3)-X(:,3).*DX(:,2)];
