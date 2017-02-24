function p = MultiPlotXY(varargin)

%MultiPlotXY - Plot two columns of each input matrix against each other.
%
% Plot the second column of each input matrix as a function of the first column.
% Optionally, alternative pairs of columns can be plotted against each other.
%
%  USAGE
%
%    p = MultiPlotXY(X,columns,Y,columns,...,<options>)
%
%    X              the first data to plot
%    columns        optional pair of columns to plot
%    Y              the second data to plot
%    columns        optional pair of columns to plot
%    ...            additional data and column indices
%    <options>      options for function <a href="matlab:help plot">plot</a>
%

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help MultiPlotXY">MultiPlotXY</a>'' for details).');
end

% We need to parse the parameters twice:
% 1) Count number of data matrices, because we need to know the number of subplots before we start plotting
n = DoIt(0,varargin{:});
% 2) Plot
DoIt(n,varargin{:});

function nData = DoIt(nSubplots,varargin)

nData = 0;
n = length(varargin);
% Parse parameters
%   i   loops through parameters
%   j   points to the first candidate parameter to be passed to 'plot'
%   k   points to the last parameter to be passed to 'plot'
i = 1;
while i <= n,
	% The first parameter is the data to plot
	X = varargin{i};
	% Default columns
	columns = [1 2];
	% Parse additional parameters if any
	if i >= n,
		j = 0;k = 0; % dummy values
	else
		% Is the second parameter a list of columns?
		if isivector(varargin{i+1},'#2','>0'),
			columns = varargin{i+1};
			j = i + 2;
		else
			j = i + 1;
		end
		% Parse additional parameters if any
		k = j-1;
		if j <= n,
			for k = j:n,
				% Is the next parameter an option for 'plot'?
				%  To answer this question, we use the following heuristics:
				%  If the next parameter is a matrix, it is *not* an option for plot, but the next data to plot,
				%  because plot does not take matrices as options - except for one possible exception: following
				%  the property 'UserData'.
				if isdmatrix(varargin{k}),
					if ~ischar(varargin{k-1}) | ~strcmp(lower(varargin{k-1}),'userdata'),
						% This parameter is a matrix and does not follow 'UserData': it is the next data to plot
						k = k - 1;
						break
					end
				end
			end
		end
	end
	nData = nData+1;
	if nSubplots ~= 0,
		subplot(nSubplots,1,nData);
		if j > i && k >= j,
			PlotXY(X,columns,varargin{j:k});
		else
			PlotXY(X,columns);
		end
	end
	i = max([i k]) + 1;
end
