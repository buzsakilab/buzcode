function [h,a] = PlotDistribution2(var1,var2,varargin)

%PlotDistribution2 - Plot bivariate data along with respective distributions.
%
%  USAGE
%
%    [h,a] = PlotDistribution2(var1,var2,<options>)
%
%    Using cell arrays will overlay variable pairs.
%
%    var1           variable 1
%    var2           variable 2
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'nBins'       number of bins for distributions (default = 100)
%     'xbins'       [m M n], lower and upper bounds, and number of bins,
%                   respectively, for x axis (default [min max 100])
%     'ybins'       [m M n], lower and upper bounds, and number of bins,
%                   respectively, for y axis (default [min max 100])
%     'smooth'      standard deviation of Gaussian kernel (default = 5)
%    =========================================================================
%
%  OUTPUT
%
%    h              handles to the plots and histograms (one line per set of
%                   variables)
%    a              axes for the plots and histograms: [main top right]


% Copyright (C) 2004-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

colors = [1 0.2 0;0 0.4 1;0 1 0.2;1 0.6 0;0.8 0 1;0.8 1 0];
% Default values
nBins = 100;
smooth = 5;
xbins = [];
ybins = [];

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotDistribution2">PlotDistribution2</a>'' for details).');
end

if isa(var1,'cell'),
	if ~isa(var2,'cell'),
		error('One parameter is a cell array and the other is not (type ''help <a href="matlab:help PlotDistribution2">PlotDistribution2</a>'' for details).');
	end
	if length(var1) ~= length(var2),
		error('Variable lists have different lengths (type ''help <a href="matlab:help PlotDistribution2">PlotDistribution2</a>'' for details).');
	end
else
	var1 = {var1};
	var2 = {var2};
end

% Check parameter sizes
for i = 1:length(var1),
	if ~isdvector(var1{i}) | ~isdvector(var2{i}) | any(size(var1{i})~=size(var2{i})),
		error('Individual variables must be vectors of identical length (type ''help <a href="matlab:help PlotDistribution2">PlotDistribution2</a>'' for details).');
	end
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotDistribution2">PlotDistribution2</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'nbins',
		nBins = varargin{i+1};
		if ~isiscalar(nBins,'>0'),
			error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help PlotDistribution2">PlotDistribution2</a>'' for details).');
		end
	case 'xbins',
		xbins = varargin{i+1};
		if ~isdvector(xbins,'#3') | xbins(1) > xbins(2) | ~isiscalar(xbins(3),'>0'),
			error('Incorrect value for property ''xbins'' (type ''help <a href="matlab:help PlotDistribution2">PlotDistribution2</a>'' for details).');
		end
	case 'ybins',
		ybins = varargin{i+1};
		if ~isdvector(ybins,'#3') | ybins(1) > ybins(2) | ~isiscalar(ybins(3),'>0'),
			error('Incorrect value for property ''ybins'' (type ''help <a href="matlab:help PlotDistribution2">PlotDistribution2</a>'' for details).');
		end
	case 'smooth',
		smooth = varargin{i+1};
		if ~isdscalar(smooth,'>=0'),
			error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help PlotDistribution2">PlotDistribution2</a>'' for details).');
		end
	otherwise,
		error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlottDistribution2">PlotDistribution2</a>'' for details).']);
  end
end

% Resize current axes to make room for distribution plots
main = gca;
position = get(main,'position');
x = position(1);
h1 = position(3); % current width
hg = h1/10; % horizontal gap
h2 = h1 - hg; % new total width
y = position(2);
v1 = position(4); % current height
vg = v1/10; % vertical gap
v2 = v1 - vg; % new total height
%  set(main,'position',[x+h2/3+hg y+v2/3+vg h2*2/3 v2*2/3]);
set(main,'position',[x y h2*2/3 v2*2/3]);
right = axes('position',[x+h2*2/3+hg y h2/3 v2*2/3]);
top = axes('position',[x y+v2*2/3+vg h2*2/3 v2/3]);

% Axis limits
if isempty(xbins),
	m1 = min(var1{1});
	M1 = max(var1{1});
	for i = 2:length(var1),
		m1 = min([m1;var1{i}]);
		M1 = max([M1;var1{i}]);
	end
	nBinsX = nBins;
else
	m1 = xbins(1);
	M1 = xbins(2);
	nBinsX = xbins(3);
end
if isempty(ybins),
	m2 = min(var2{1});
	M2 = max(var2{1});
	for i = 2:length(var1),
		m2 = min([m2;var2{i}]);
		M2 = max([M2;var2{i}]);
	end
	nBinsY = nBins;
else
	m2 = ybins(1);
	M2 = ybins(2);
	nBinsY = ybins(3);
end
d1 = M1-m1;
m1 = m1-d1/10;
M1 = M1+d1/10;
d2 = M2-m2;
m2 = m2-d2/10;
M2 = M2+d2/10;

h = [];
a = [main top right];

for i = 1:length(var1),
	% Plot bivariate data
	axes(main);hold on;
	p = plot(var1{i},var2{i},'o');
	k = mod(i,size(colors,1)+1);
	set(p,'MarkerFaceColor',colors(k,:),'MarkerEdgeColor',colors(i,:));
	if i > 1,
		m1 = min([m1 pm1]);
		M1 = max([M1 pM1]);
		m2 = min([m2 pm2]);
		M2 = max([M2 pM2]);
	end
	pm1 = m1;
	pM1 = M1;
	pm2 = m2;
	pM2 = M2;
	xlim([m1 M1]);ylim([m2 M2]);

	% Distribution of var1
	axes(top);hold on;
	[dist1,x1] = hist(var1{i},m1:(M1-m1)/nBinsX:M1);
	dist1 = Smooth(dist1/sum(dist1),smooth);
	bv = bar(x1,dist1,1);
	k = mod(i,size(colors,1)+1);
	set(bv,'FaceColor',colors(k,:),'EdgeColor',colors(i,:));
	set(top,'xtick',[]);
	xlim([m1 M1]);

	% Distribution of var2
	axes(right);hold on;
	[dist2,x2] = hist(var2{i},m2:(M2-m2)/nBinsY:M2);
	dist2 = Smooth(dist2/sum(dist2),smooth);
	bh = barh(x2,dist2,1);
	k = mod(i,size(colors,1)+1);
	set(bh,'FaceColor',colors(k,:),'EdgeColor',colors(i,:));
	set(right,'ytick',[]);
	ylim([m2 M2]);
	
	% store
	h = [h;p bv bh];
end

% Adjust axis locations
set(main,'xaxislocation','bottom','yaxislocation','left');
set(right,'xaxislocation','bottom','yaxislocation','right');
set(top,'xaxislocation','top','yaxislocation','left');
axes(main);