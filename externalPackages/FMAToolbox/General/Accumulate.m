function [a,b,c] = Accumulate(variables,values,outputSize)

%Accumulate - Accumulate repeated observations.
%
%  Accumulate repeated observations Y(i) = f(X1(i),X2(i),...,Xn(i)), where
%  the Xj(i) are binned variables (i.e. they run from 1 to Nj). This
%  is useful for instance to compute:
%
%   * a PSTH, where Y(i) = 1 and X(i) lists the binned spike timestamps
%     relative to the synchronizing events
%   * a JPSTH, where Y(i) = 1, X1(i) lists the binned spike timestamps of
%     the first neuron, and X2(i) lists those of the second neuron
%   * the dwell time of the animal in the environment, where Y(i) is the
%     sampling interval, and X1(i) and X2(i) list the binned x and y
%     coordinates of the animal for each video sample
%   * the spatial spike count of a place cell, where Y(i) = 1, and X1(i) and
%     X2(i) list the binned x and y coordinates of the animal for each spike
%   * the firing map of a place cell, as the ratio of the spatial spike count
%     and the dwell time
%   * the angular firing curve of a head direction cell, as the ratio of the
%     angular spike count and the dwell time
%   * etc.
%
%  USAGE
%
%    % Sum repeated observations
%    sum = Accumulate(variables,values,outputSize)
%
%    % Average repeated observations, and compute std and 95% c.i.
%    % NOTE: unbiased std and confidence intervals are valid only for linear data
%    [mean,std,conf] = Accumulate(variables,values,outputSize)
%
%    variables      M x N matrix, listing the M instances of the N variables
%                   X1(i)...Xn(i)
%    values         optional M instances of Y(i) (circular values (e.g. angles)
%                   must be provided in complex exponential form exp(i*theta),
%                   where theta is in radians) (default = 1)
%    outputSize     optional output matrix size (default = max(variables))
%
%  EXAMPLE
%
%    X = Bin(x,[0 1],nBins);       % bin x coordinate
%    Y = Bin(y,[0 1],nBins);       % bin y coordinate
%    time = Accumulate([X Y],dt);  % compute occupancy map
%
%    (This is only a simplified example; to actually compute an occupancy map,
%    use FiringMap instead)
%
%  NOTES
%
%    1) Angular data should be provided in complex exponential form
%    2) Standard error and confidence intervals are valid only for linear data
%
%  SEE
%
%    See also Bin.

% Copyright (C) 2004-2011 by Michaël Zugaro, 2004 by Ken Harris
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

a = [];
b = [];
c = [];

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Accumulate">Accumulate</a>'' for details).');
end

if isempty(variables),
	return;
end

if nargin < 2
	values = 1;
end

if nargin < 3
	outputSize = max(variables);
end

% Make sure 'values' is a Mx1 vector (possibly a scalar)
if ~isdvector(values),
    error('Incorrect values (type ''help <a href="matlab:help Accumulate">Accumulate</a>'' for details)');
end
if size(values,1) == 1, values = values'; end

% Make sure 'variables' is a MxN matrix (possibly a Mx1 vector) of positive integers
if ~isivector(variables,'>=0') & ~isimatrix(variables,'>=0'),
    error('Incorrect variables (type ''help <a href="matlab:help Accumulate">Accumulate</a>'' for details)');
end
if size(variables,1) == 1, variables = variables'; end

% If 'values' is a scalar, make it a vector the same length as 'variables'
if length(values) == 1,
	values = repmat(values,size(variables,1),1);
end

% If 'variables' is a vector, make it a matrix
if size(variables,2) == 1,
	variables = [variables(:) ones(length(variables),1)];
	if length(outputSize) == 1,
		outputSize = [outputSize 1];
	end
end

% 'variables' and 'values' must have the same length
if size(values,1) ~= size(variables,1),
    error('Incompatible sizes for variables and values (type ''help <a href="matlab:help Accumulate">Accumulate</a>'' for details)');
end

% Drop NaN and Inf values
i = any(isnan(variables)|isinf(variables),2)|isnan(values)|isinf(values);
variables(i,:) = [];
values(i,:) = [];

% No variable should be outside the output size range
if any(max(variables,[],1)>outputSize),
	error('Variables are out of bounds (type ''help <a href="matlab:help Accumulate">Accumulate</a>'' for details)');
end

% Check output size
if ~isdvector(outputSize,'>0') || length(outputSize) ~= size(variables,2),
	error('Incorrect output size (type ''help <a href="matlab:help Accumulate">Accumulate</a>'' for details)');
end

% Make linear index from subscripts
for i = 1:size(variables,2),
	subscript{i} = variables(:,i);
end
linearIndex = sub2ind(outputSize,subscript{:});

tmp = sparse(linearIndex,1,values,prod(outputSize),1);
a = reshape(full(tmp),outputSize);

if nargout >= 2,
	% Mean: divide by cell count
	tmp = sparse(linearIndex,1,1,prod(outputSize),1);
	n = reshape(full(tmp),outputSize);
	a = a ./ n;
	a(isnan(a)) = 0;
	% Standard deviation: use V = E(X2)-E(X)2
	tmp = sparse(linearIndex,1,values.^2,prod(outputSize),1);
	s = reshape(full(tmp),outputSize)./n; % E(X2)
	b = sqrt((s-a.^2)./(n-1).*n); % use the unbiased formula
	b(isnan(b)) = 0;
	% 95% confidence intervals
	if length(size(variables)) <= 2,
		% not computed for > 2 dimensions (too long)
		if isvector(a),
			for i = 1:length(a),
				ok = variables(:,1) == i;
				c(i,1) = prctile(values(ok),5);
				c(i,2) = prctile(values(ok),95);
			end
		else
			for i = 1:size(a,1),
				for j = 1:size(a,2),
					ok = variables(:,1) == i & variables(:,2) == j;
					c(i,j,1) = prctile(values(ok),5);
					c(i,j,2) = prctile(values(ok),95);
				end
			end
		end
	end
end
