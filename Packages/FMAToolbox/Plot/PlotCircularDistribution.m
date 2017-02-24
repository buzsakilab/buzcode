function PlotCircularDistribution(dist,a,s,varargin)

%PlotCircularDistribution - Plot phase distribution and statistics.
%
%  Plot one or more phase distributions as histograms and polar plots
%  showing circular mean and resultant length. Dotted vectors indicate
%  that the data is not significantly clustered around a single peak
%  (Rayleigh test, alpha > 0.05).
%
%  USAGE
%
%    PlotCircularDistribution(dist,angles,stats)
%
%    dist           phase distribution(s) (obtained using <a href="matlab:help CircularDistribution">CircularDistribution</a>)
%    angles         optional centers of the angular bins
%    stats          optional distribution statistics
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'mode'        one plot for all distributions ('grouped', default) or
%                   one per distribution ('separate')
%    =========================================================================
%
%  SEE
%
%    See also CircularDistribution.

% Copyright (C) 2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

mode = 'grouped';
v = varargin;

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotCircularDistribution">PlotCircularDistribution</a>'' for details).');
end

% Check parameter size
if ~isdvector(dist) && ~isdmatrix(dist),
	error('Incorrect distributions (type ''help <a href="matlab:help PlotCircularDistribution">PlotCircularDistribution</a>'' for details).');
end

% Stats provided?
if nargin >= 3 && ~isstruct(s),
	v = {s,varargin{:}};
end
if nargin < 3 || ~isstruct(s),
	stats = [];
else
	stats = s;
end

% Angles provided?
if nargin >= 2 && ischar(a),
	v = {a,s,varargin{:}};
end
if nargin < 2 || ischar(a),
	angles = linspace(0,2*pi,size(dist,1)+1)';angles(end) = [];
	binSize = angles(2)-angles(1);
	angles = angles + binSize/2;
else
	angles = a;
end

varargin = v;

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotCircularDistribution">PlotCircularDistribution</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'mode',
			mode = lower(varargin{i+1});
			if ~isstring_FMAT(mode,'separate','grouped'),
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help PlotCircularDistribution">PlotCircularDistribution</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotCircularDistribution">PlotCircularDistribution</a>'' for details).']);
	end
end

nGroups = size(dist,2);
color = Bright(nGroups);

% Plot histograms (or curves) for each group
M = 0; % max y axis value
for i = 1:nGroups,
	h = dist(:,i);
	if nGroups == 1,
		subplot(1,2,1); hold on;
		bar([angles;angles+2*pi],[h;h],1,'EdgeColor','none','FaceColor','b');
	else
		if strcmp(mode,'grouped'),
			subplot(1,2,1);
		else
			subplot(2,nGroups,i);
		end
		hold on;
		plot([angles;angles+2*pi],[h;h],'color',color(i,:));
	end
	m = ylim;m = m(2);
	M = max([m M]);
end
% Adjust axes and labels, plot sinewave
if strcmp(mode,'separate'),
	m = 2;
	n = nGroups;
	N = nGroups;
else
	m = 1;
	n = 2;
	N = 1;
end
for i = 1:N,
	subplot(m,n,i);
	xlim([0 4*pi]);
	x = linspace(0,4*pi,600);
	a = M/4;
	y = a*cos(x)+M+a;
	plot(x,y,'r','linewidth',2);
	xlabel('Angle (rad)');
	ylabel('Probability');
end

% Show statistics for each group
for i = 1:nGroups,
	if strcmp(mode,'grouped'),
		subplot(1,2,2);
	else
		subplot(2,nGroups,i+nGroups);
	end
	hold on;
	h = dist(:,i);
	p = polar([angles;angles(1)],[h;h(1)]);
	set(p,'color',color(i,:));
	if ~isempty(stats),
		c = compass(stats.r(i)*cos(stats.m(i)),stats.r(i)*sin(stats.m(i)));
		set(c,'color',color(i,:));
		if stats.p(i) > 0.05, set(c,'linestyle',':'); end
	end
end
for i = 1:N,
	subplot(m,n,i+N);
	axis square;
	plot([-M M],[0 0],'k:');
	plot([0 0],[-M M],'k:');
	xlim([-M M]);ylim([-M M]);
	set(gca,'xtick',[],'ytick',[],'box','on');
	if i == 1,
		xlabel('Circular Mean and Resultant Length');
	end
end