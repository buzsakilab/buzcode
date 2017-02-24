function PlotShortTimeCCG(data,varargin)

%PlotShortTimeCCG - Plot time-varying auto/cross-correlograms of point processes.
%
%  USAGE
%
%    PlotShortTimeCCG(data,<options1>,<options2>)
%
%    data           data obtained using <a href="matlab:help ShortTimeCCG">ShortTimeCCG</a>
%    <options1>     optional list of property-value pairs (see table below)
%    <options2>     optional list of property-value pairs (see <a href="matlab:help PlotColorMap">PlotColorMap</a>)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'parent'      parent figure or uipanel handle (default = gcf)
%    =========================================================================
%
%  SEE
%
%    See also CCG, ShortTimeCCG.

% Copyright (C) 2010-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
parent = [];
x = 1:size(data,2);
y = 1:size(data,1);
args = {};

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help PlotShortTimeCCG">PlotShortTimeCCG</a>'' for details).');
end

% Parse parameter list
j = 1;
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotShortTimeCCG">PlotShortTimeCCG</a>'' for details).']);
	end
	switch(lower(varargin{i})),

		case 'parent',
			parent = varargin{i+1};
			if ~ishandle(parent),
				error('Incorrect value for property ''parent'' (type ''help <a href="matlab:help PlotShortTimeCCG">PlotShortTimeCCG</a>'' for details).');
			end

		case 'x',
			x = varargin{i+1};
			if ~isdvector(x),
				error('Incorrect value for property ''x'' (type ''help <a href="matlab:help PlotShortTimeCCG">PlotShortTimeCCG</a>'' for details).');
			end
			args{j} = varargin{i};
			args{j+1} = varargin{i+1};
			j = j + 2;

		case 'y',
			y = varargin{i+1};
			if ~isdvector(y),
				error('Incorrect value for property ''y'' (type ''help <a href="matlab:help PlotShortTimeCCG">PlotShortTimeCCG</a>'' for details).');
			end
			args{j} = varargin{i};
			args{j+1} = varargin{i+1};
			j = j + 2;

		case {'threshold','cutoffs','gamma','colorspace','bar','type'},
			args{j} = varargin{i};
			args{j+1} = varargin{i+1};
			j = j + 2;

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotShortTimeCCG">PlotShortTimeCCG</a>'' for details).']);

	end
end

if isempty(parent), parent = gcf; end

xcorrs = round(linspace(1,size(data,2),5));
for i = 1:5,
	subplot(5,6,6*i,'parent',parent);
	j = round(xcorrs(i));
	bar(y,data(:,j));
	xlim([min(y) max(y)]);
	xlabel(['t = ' num2str(x(j))]);
	set(gca,'ytick',[]);
end
subplot(5,6,[1 29],'parent',parent);
PlotColorMap(data,1,args{:});
xlabel('Time (s)');
ylabel('Time shift (s)');

