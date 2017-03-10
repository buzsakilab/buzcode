function [bins,binned] = Bin(x,a,b,c)

%Bin - Assign each input value a bin number between 1 and N.
%
%  Assign each input value a bin number between 1 and N. Also return the
%  corresponding bin values (lower bounds), e.g. assuming that the bins
%  are [0,1] [1,2] [2,3] ... the value 2.3 falls in bin #3 and the bin
%  value (lower bound) is 2.
%
%  USAGE
%
%    % To specify bin limits and # bins
%    [bins,binned] = Bin(x,limits,nBins)
%
%    % To automatically set the limits to [min(x) max(x)]
%    bins = Bin(x,nBins)
%
%    % To specify an explicit list of bins (including the upper bound
%    % of the last bin)
%    [bins,binned] = Bin(x,bins)
%
%    % To discard values out of limits
%    bins = Bin(x,...,'trim')
%
%    x              1D variable
%    limits         [min_x max_x], where min_x is included and max_x is excluded;
%                   these need not be equal to min(x) and max(x)
%                   by default, x values out of bounds are counted in the first
%                   and last bin, respectively
%    nBins          number of bins
%    bins           explicit list of bins - including the upper bound of the
%                   last bin
%
%  NOTE
%
%    y = Bin(x,[m M],N) uses N bins b(i). Their size is k=(M-m)/N and their
%    lower and upper bounds are m+i*k (included) and m+(i+1)*k (excluded),
%    respectively. As intended, the extreme bounds are equal to m and M,
%    respectively, but keep in mind that while the first bin starts at m,
%    the last bin *ends* at M (it actually starts at m+(N-1)*k).
%
%  EXAMPLES
%
%    Bin([2 2.3 4 4.7],[2 5],3);  % returns [1 1 3 3]
%    Bin([2 2.3 4 4.7],[2 5],6);  % returns [1 1 5 6]
%    Bin([2 2.3 4 4.2 5],[2 5],3);  % returns [1 1 3 3 3]
%
%  SEE
%
%    See also Accumulate.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Bin">Bin</a>'' for details).');
end

if isempty(x),
	bins = [];
	binned = [];
	return
end
if ~isdvector(x),
	error('Data is not a vector (type ''help <a href="matlab:help Bin">Bin</a>'' for details).');
end

%  if ~isdvector(a),
%  	error('Incorrect second parameter (type ''help <a href="matlab:help Bin">Bin</a>'' for details).');
%  end

if length(a) == 1,
	% Second parameter is a scalar, i.e. the number of bins
	if ~isiscalar(a,'>0'),
		error('Incorrect number of bins (type ''help <a href="matlab:help Bin">Bin</a>'' for details).');
	end
	nBins = a;
	% Automatic limits
	limits = [min(x) max(x)];
	% Trim?
	trim = nargin == 3 && strcmp(lower(b),'trim');
elseif length(a) == 2,
	% Second parameter is a pair, i.e. the limits
	limits = a;
	% Third parameter must be # bins
	if ~isiscalar(b,'>0'),
		error('Incorrect number of bins (type ''help <a href="matlab:help Bin">Bin</a>'' for details).');
	end
	nBins = b;
	% Trim?
	trim = nargin == 4 && strcmp(lower(c),'trim');
else
	% Explicit bins
	d = a(2) - a(1);
%  	if any(abs(diff(a)-d)>eps),
%  		error('Incorrect list of bins - bins must be evenly spaced (type ''help <a href="matlab:help Bin">Bin</a>'' for details).');
%  	end
	% Determine limits and # bins
	limits = [a(1) a(end)];
	nBins = length(a)-1;
	% Trim?
	trim = nargin == 3 && strcmp(lower(b),'trim');
end

if trim,
	x(x<limits(1)) = nan;
	x(x>=limits(2)) = nan;
else
	x(x<limits(1)) = limits(1);
	x(x>=limits(2)) = limits(2);
end

bins = floor((x-limits(1))/(limits(2)-limits(1))*nBins)+1;

if ~trim,
	bins(bins>nBins) = nBins;
end

binSize = (limits(2)-limits(1))/nBins;
i = linspace(limits(1),limits(2),nBins+1);
n = isnan(bins);
binned(~n,1) = i(bins(~n))+binSize/2;
binned(n,1) = NaN;
