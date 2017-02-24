function PlotPhaseDistribution(dist,angles,stats)

%PlotPhaseDistribution - Plot phase distribution.
%
%  Plot one or more phase distributions as histograms and polar plots
%  showing circular mean and resultant length. Dotted vectors indicate
%  that the data is not significantly clustered around a single peak
%  (Rayleigh test, alpha > 0.05).
%
%  USAGE
%
%    PlotPhaseDistribution(dist,angles,stats)
%
%    dist           phase distribution(s) (obtained using <a href="matlab:help PhaseDistribution">PhaseDistribution</a>)
%    angles         optional centers of the angular bins
%    stats          optional distribution statistics
%
%  SEE
%
%    See also PhaseDistribution.

% Copyright (C) 2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotPhaseDistribution">PlotPhaseDistribution</a>'' for details).');
end

% Check parameter size
if ~isdvector(dist) && ~isdmatrix(dist),
	error('Incorrect distributions (type ''help <a href="matlab:help PlotPhaseDistribution">PlotPhaseDistribution</a>'' for details).');
end

% Angles
if nargin < 2,
	angles = linspace(0,2*pi,size(dist,1)+1);angles(end) = [];
	binSize = angles(2)-angles(1);
	angles = angles + binSize/2;
end

% Plot histograms (or curves) for each group
if nargin == 3,
	subplot(1,2,1);
end
hold on;
nGroups = size(dist,2);
color = Bright(nGroups);
for i = 1:nGroups,
	h = dist(:,i);
	if nGroups == 1,
		bar([angles;angles+2*pi],[h;h],1,'EdgeColor','none','FaceColor','b');
	else
		plot([angles;angles+2*pi],[h;h],'color',color(i,:));
	end
end

% Adjust axes and labels, plot sinewave
xlim([0 4*pi]);
M = ylim;M = M(2);
x = linspace(0,4*pi,600);
a = M/4;
y = a*cos(x)+M+a;
plot(x,y,'r','linewidth',2);
xlabel('Angle (rad)');
ylabel('Probability');

% Show statistics for each group
if nargin == 3,
	subplot(1,2,2);
	hold on;
	for i = 1:nGroups,
		h = dist(:,i);
		p = polar([angles;angles(1)],[h;h(1)]);
		set(p,'color',color(i,:));
		%p = polar([0 stats.m(i)],[0 stats.r(i)]);
		%set(p,'color',color(i,:));
		c = compass(stats.r(i)*cos(stats.m(i)),stats.r(i)*sin(stats.m(i)));
		set(c,'color',color(i,:));
		if stats.p(i) > 0.05, set(c,'linestyle',':'); end
	end
	axis square;
	plot([-M M],[0 0],'k:');
	plot([0 0],[-M M],'k:');
	xlim([-M M]);ylim([-M M]);
	set(gca,'xtick',[],'ytick',[],'box','on');
	xlabel('Circular Mean and Resultant Length');
end
