function p = PlotMean(X,Y,L,U,style,color)

%PlotMean - Plot mean and confidence intervals.
%
%  USAGE
%
%    p = PlotMean(X,Y,L,U,style,color)
%
%    X              abscissae
%    Y              ordinates
%    L              lower error level
%    U              upper error level
%    style          optional style (':' or '-')
%    color          optional color specification (see <a href="matlab:help plot">plot</a>)
%
%  NOTES
%
%    Notice that L and U specifiy actual error levels rather than (relative)
%    ranges, contrary to the built-in 'errorbar' where the error is in [Y-L Y+U].
%
%    If X is empty, it is set to 1:N (where N is the length of Y). If Y is
%    empty, only error limits are plotted. Conversely, if L and U are empty, only
%    the mean is plotted.

% Copyright (C) 2008-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 4,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotMean">PlotMean</a>'' for details).');
end

if nargin < 6,
	color = 'b';
end
if nargin < 5,
	style = ':';
end

nX = length(X);
nY = length(Y);
nL = length(L);
nU = length(U);

% Make sure all non-empty inputs have the same length
n = [nX nY nL nU];
n = n(n~=0);
if isempty(n),
	error('At least one input must contain data (type ''help <a href="matlab:help PlotMean">PlotMean</a>'' for details).');
end
m = n - n(1);
if any(m),
	error('All non-empty inputs must have the same length (type ''help <a href="matlab:help PlotMean">PlotMean</a>'' for details).');
end

if nX == 0,
	X = (1:n(1))';
end

hold on;
% Plot mean
if nY ~= 0,
	p = plot(X,Y);
	if strcmp(style,':'),
		set(p,'LineStyle','-','color',color);
	else
		set(p,'LineWidth',3,'color',color);
	end
end
% Plot lower error level
if nL ~= 0,
	p = plot(X,L);
	if strcmp(style,':'),
		set(p,'LineStyle',':','color',color);
	else
		set(p,'LineWidth',1,'color',color);
	end
end
% Plot upper error level
if nU ~= 0,
	p = plot(X,U);
	if strcmp(style,':'),
		set(p,'LineStyle',':','color',color);
	else
		set(p,'LineWidth',1,'color',color);
	end
end

