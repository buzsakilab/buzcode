function p = PlotXY(X,varargin)

%PlotXY - Plot two columns of a matrix against each other.
%
% Plot the second column of a matrix as a function of the first column.
% Optionally, alternative pairs of columns can be plotted against each other.
%
%  USAGE
%
%    p = PlotXY(X,columns,<options>)
%
%    X              the data to plot
%    columns        optional pair of columns to plot
%    <options>      options for function <a href="matlab:help plot">plot</a>
%

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotXY">PlotXY</a>'' for details).');
end
if nargin >= 2 & isdvector(varargin{1},'#2','>0'),
	% Parameter #2 is a pair of column numbers
	if size(X,2) == 1,
		p = plot(X,varargin{2:end});
	else
		p = plot(X(:,varargin{1}(1)),X(:,varargin{1}(2)),varargin{2:end});
	end
else
	if size(X,2) == 1,
		p = plot(X,varargin{:});
	else
		p = plot(X(:,1),X(:,2),varargin{:});
	end
end
