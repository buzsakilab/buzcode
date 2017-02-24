function p = npcdf(X,x)

%npcdf - Non-parametric cumulative distribution function.
%
% Given a list of samples drawn from an unknown (non-parametric) distribution,
% the corresponding cumulative distribution function is evaluated at a set of
% given points.
%
%  USAGE
%
%    p = npcdf(X,x)
%
%    X              list of samples drawn from the unknown distribution
%    x              set of points where the cdf should be evaluated
%

% Copyright (C) 2010-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help npcdf">npcdf</a>'' for details).');
end

% Get rid of NaNs, reshape inputs as vectors
X(isnan(X)) = [];
X = X(:);
nX = length(X);
x(isnan(x)) = [];
x = x(:);
nx = length(x);

% Construct a matrix where the first column is a concatenation of the samples and points,
% the second is 1 for samples and 0 for points, and the third is the order of the points
% (this will be necessary to keep track of the points after the matrix is reordered)
Y = [X ones(nX,1) zeros(nX,1);x zeros(nx,1) (1:nx)'];
% Sort samples and points in ascending order, so that the position of the points in the matrix
% (row numbers) corresponds to their cdf values
[Y,i] = sortrows(Y);

% Actually, row numbers should only take samples into account; this is precisely what the
% second column of Y is for (it has 1 for samples and 0 for points)
F = cumsum(Y(:,2))/nX;
% Now that we have the value of the cdf at each sample and each point, we need to extract
% the value at the points; again, we use the second column of Y (it has 0 for points)
tested = Y(:,2) == 0;
% Because Y was reordered (sortrows), we need to reorder the values of F back to the order
% in which the points were listed in x
p(Y(tested,3)) = F(tested);
