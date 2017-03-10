function PlotPhasePrecession(data,stats,varargin)

%PlotPhasePrecession - Plot phase precession plots and maps.
%
%  USAGE
%
%    PlotPhasePrecession(data,stats,<options>)
%
%    data           phase precession data obtained using <a href="matlab:help PhasePrecession">PhasePrecession</a>
%    stats          optional phase precession stats obtained using <a href="matlab:help PhasePrecession">PhasePrecession</a>
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      size in bins of the smoothing kernel for the phase map
%                   (0 = no smoothing, default = 4)
%     'nBins'       number of x and phase bins (default = [100 100])
%     'track'       'linear' or 'circular' (default = 'linear')
%     'laps'        plot all laps together ('all', default) or separately
%                   ('single')
%     'parent'      parent figure or uipanel handle (default = gcf)
%    =========================================================================
%
%  SEE
%
%    See also PhasePrecession.

% Copyright (C) 2009-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
track = 'linear';
nBins = [100 100];
smooth = 4;
parent = [];
laps = 'all';

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).');
end

if nargin > 1 && ~isstruct(stats),
	varargin = {stats varargin{:}};
	stats = [];
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'smooth',
			smooth = varargin{i+1};
			if ~isdscalar(smooth,'>0'),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).');
			end
		case 'nbins',
			nBins = varargin{i+1};
			if ~isivector(nBins,'>0','#2'),
				error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).');
			end
		case 'track',
			track = lower(varargin{i+1});
			if ~isstring_FMAT(track,'linear','circular'),
				error('Incorrect value for property ''track'' (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).');
			end
		case 'laps',
			laps = lower(varargin{i+1});
			if ~isstring_FMAT(laps,'all','single'),
				error('Incorrect value for property ''laps'' (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).');
			end
		case 'parent',
			parent = varargin{i+1};
			if ~ishandle(parent),
				error('Incorrect value for property ''parent'' (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).']);
	end
end

if isempty(parent), parent = gcf; end

% Initialize a few variables
t = data.position.t;
x = data.position.x;
pp = data.position.phase;
l = data.position.lap;
r = data.rate.r;
pr = data.rate.phase;
nBinsX = nBins(1);
nBinsP = nBins(2);

if strcmp(laps,'all'),

	% Rate vs position plot
	a = subplot(4,1,1,'parent',parent);
	nBins = 200;
	curve = FiringCurve(data.x,t,'nbins',nBins,'smooth',2,'type','c');
	if strcmp(track,'linear'),
		plot(curve.x,curve.rate,'b');
		xlim([0 1]);
	else
		plot([curve.x curve.x+1],[curve.rate curve.rate],'b');
		xlim([0 2]);
	end
	set(a,'tickdir','out');
	ylabel('Firing Rate (Hz)');

	% Phase vs position plot
	a = subplot(4,1,2,'parent',parent);
	if strcmp(track,'linear'),
		plot([x;x],[pp;pp+2*pi]*180/pi,'k.','markersize',4);
		if ~isempty(stats),
			hold on;
			dx = max(x)-min(x);
			xm = mean([min(x);max(x)]);
			ym = (stats.slope*xm+stats.intercept)*180/pi;
			PlotSlope(xm,ym,stats.slope*180/pi,dx,'r');
			PlotSlope(xm,ym+360,stats.slope*180/pi,dx,'r');
		end
		xlim([0 1]);
	else
		plot([x;x;x+1;x+1],[pp;pp+2*pi;pp;pp+2*pi]*180/pi,'k.','markersize',4);
		if ~isempty(stats),
			hold on;
			dx = max(x)-min(x);
			xm = mean([min(x);max(x)]);
			ym = (stats.slope*xm+stats.intercept)*180/pi;
			PlotSlope(xm,ym,stats.slope*180/pi,dx,'r');
			PlotSlope(xm,ym+360,stats.slope*180/pi,dx,'r');
			PlotSlope(xm+1,ym,stats.slope*180/pi,dx,'r');
			PlotSlope(xm+1,ym+360,stats.slope*180/pi,dx,'r');
		end
		xlim([0 2]);
	end
	set(a,'ytick',0:90:720,'ylim',[0 720],'tickdir','out');
	ylabel('Phase (°)');

	% Phase vs position map
	a = subplot(4,1,3,'parent',parent);
	x0 = Bin(x,[0 1],nBinsX);
	p0 = Bin(pp,[0 2*pi],nBinsP);
	map = Accumulate([p0 x0],1,[nBinsP nBinsX]);
	if strcmp(track,'linear'),
		map = [map;map];
		map = Smooth(map,2);
		PlotColorMap(map,1,'x',linspace(0,1,nBinsX+1),'y',0:6:720);
	else
		map = [map map;map map];
		map = Smooth(map,smooth);
		PlotColorMap(map,1,'x',linspace(0,2,2*nBinsX+1),'y',0:6:720);
	end
	set(a,'ytick',0:90:720,'ylim',[0 720],'tickdir','out','ydir','normal');
	xlabel('Position');
	ylabel('Phase (°)');

	% Phase vs rate plot
	a = subplot(4,1,4,'parent',parent);
	if isempty(r), return; end
	rnd = randn(size(r))/5; % jitter points horizontally
	plot(r+rnd,pr*180/pi,'k.','markersize',4);
	hold on;
	n = max(r);
	for i = 1:n,
		if any(r==i),
			[m(i,1),c(i,:)] = CircularConfidenceIntervals(pr(r==i));
		else
			m(i,1) = NaN;
			c(i,1:2) = NaN;
		end
	end
	neg = m<0;
	m(neg) = m(neg)+2*pi;
	c(neg,:) = c(neg,:)+2*pi;
	m = m*180/pi;
	c = c*180/pi;
	errorbar(1:n,m,m-c(:,2),c(:,1)-m,'r');
	set(a,'ytick',0:90:360,'ylim',[0 360],'tickdir','out');
	xlabel('Firing Rate (spikes / 2 cycles)');
	ylabel('Phase (°)');

else

	colors = 'bk';

	% List and count lap numbers
	laps = unique(data.position.lap);
	laps = laps(laps~=0);
	nLaps = length(laps);

	% Rate vs position plot
	a = subplot(2,1,1,'parent',parent);
	hold on;
	nBins = 200;
	curve = FiringCurve(data.x,t,'nbins',nBins,'smooth',2,'type','c');
	for i = 1:nLaps,
		style = colors(mod(i,2)+1);
		if strcmp(track,'linear'),
			plot(curve.x+(i-1),curve.rate,style,'markersize',4);
			xlim([0 nLaps]);
		else
			plot([curve.x curve.x+1]+2*(i-1),[curve.rate curve.rate],style,'markersize',4);
			xlim([0 2*nLaps]);
		end
	end
	set(a,'tickdir','out','xtick',1:nLaps,'xticklabel',[]);
	xlabel('Position');
	ylabel('Firing Rate (Hz)');

	% Phase vs position plot
	a = subplot(2,1,2,'parent',parent);
	hold on;
	for i = 1:nLaps,
		style = [colors(mod(i,2)+1) '.'];
		if strcmp(track,'linear'),
			plot([x(l==i);x(l==i)]+(i-1),[pp(l==i);pp(l==i)+2*pi]*180/pi,style,'markersize',4);
			if ~isempty(stats),
				hold on;
				dx = max(x(l==i))-min(x(l==i));
				xm = mean([min(x(l==i));max(x(l==i))])+(i-1);
				ym = (stats.lap.slope(i)*xm+stats.lap.intercept(i))*180/pi;
				PlotSlope(xm,ym,stats.lap.slope(i)*180/pi,dx,'r');
				PlotSlope(xm,ym+360,stats.lap.slope(i)*180/pi,dx,'r');
			end
			xlim([0 nLaps]);
		else
			plot([x(l==i);x(l==i);x(l==i)+1;x(l==i)+1]+2*(i-1),[pp(l==i);pp(l==i)+2*pi;pp(l==i);pp(l==i)+2*pi]*180/pi,style,'markersize',4);
			if ~isempty(stats),
				hold on;
				dx = max(x(l==i))-min(x(l==i));
				xm = mean([min(x(l==i));max(x(l==i))]);
				ym = (stats.lap.slope(i)*xm+stats.lap.intercept(i))*180/pi;
				PlotSlope(xm+2*(i-1),ym,stats.lap.slope(i)*180/pi,dx,'r');
				PlotSlope(xm+2*(i-1),ym+360,stats.lap.slope(i)*180/pi,dx,'r');
				PlotSlope(xm+1+2*(i-1),ym,stats.lap.slope(i)*180/pi,dx,'r');
				PlotSlope(xm+1+2*(i-1),ym+360,stats.lap.slope(i)*180/pi,dx,'r');
			end
			xlim([0 2*nLaps]);
		end
	end
	set(a,'ytick',0:90:720,'ylim',[0 720],'tickdir','out','xtick',1:nLaps,'xticklabel',[]);
	xlabel('Position');
	ylabel('Phase (°)');

end