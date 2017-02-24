function PlotCSD(csd,varargin)

%PlotCSD - Plot current source density.
%
%  USAGE
%
%    PlotCSD(csd,<options>)
%
%    csd            current source density (see <a href="matlab:help CSD">CSD</a>)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'lfp'         local field potential
%     'scale'       scale (arbitrary units) for LFP traces (default = 1)
%     'cutoffs'     cutoff values (default = [-M M] where M is the maximum
%                   amplitude of the CSD)
%    =========================================================================
%
%  SEE
%
%    See <a href="matlab:help CSD">CSD</a>.

% Copyright (C) 2008-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
lfp = [];
scale = 1;
cutoffs = [];

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotCSD">PlotCSD</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help PlotCSD">PlotCSD</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'lfp',
			lfp = varargin{i+1};
			if ~isdmatrix(lfp) | size(lfp,1) ~= size(csd,1) | size(lfp,2) ~= size(csd,2)+2,
				error('Incorrect size for parameter ''lfp'' (type ''help <a href="matlab:help PlotCSD">PlotCSD</a>'' for details).');
			end
		case 'scale',
			scale = varargin{i+1};
			if ~isdscalar(scale,'>0'),
				error('Incorrect value for property ''scale'' (type ''help <a href="matlab:help PlotCSD">PlotCSD</a>'' for details).');
			end
		case 'cutoffs',
			cutoffs = varargin{i+1};
			if ~isdvector(cutoffs,'#2','<'),
				error('Incorrect value for property ''cutoffs'' (type ''help <a href="matlab:help PlotCSD">PlotCSD</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotCSD">PlotCSD</a>'' for details).']);
	end
end

d = csd(:,2:end);
d = interp2(d);
d = d(1:2:size(d,1),:);

pcolor(csd(:,1),1:size(d,2),flipud(transpose(d)));
shading interp;
if ~isempty(cutoffs),
	set(gca,'clim',cutoffs);
end

if ~isempty(lfp),
	hold on;
	y = lfp(:,2:end);
	y = y - repmat(mean(y),size(y,1),1);
	y = y / max(max(abs(y)))*scale;
	n = size(y,2);
	for i = 1:n,
		plot(lfp(:,1),y(:,i)+(n-i)*2-1,'k');
    end
    ylim([-2 (n-1)*2]);
end

