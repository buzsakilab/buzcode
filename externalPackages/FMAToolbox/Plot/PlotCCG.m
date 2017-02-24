function PlotCCG(t,ccg,varargin)

%PlotCCG - Plot auto/cross-correlograms of point processes.
%
%  USAGE
%
%    PlotCCG(t,ccg)
%
%    t              time bins obtained using <a href="matlab:help CCG">CCG</a>
%    ccg            data obtained using <a href="matlab:help CCG">CCG</a>
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'full'        'on' to plot all correlograms, 'off' to plot only upper
%                   triangular matrix of correlograms (default)
%    =========================================================================
%
%  SEE
%
%    See also CCG, CCGParameters, ShortTimeCCG.

% Copyright (C) 2010-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
full = false;

if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help PlotCCG">PlotCCG</a>'' for details).');
end

if ~isdvector(t),
	error('Incorrect times bins (type ''help <a href="matlab:help PlotCCG">PlotCCG</a>'' for details).');
end
if length(t) ~= size(ccg,1),
	error('Inconsistent time bins and cross-correlogram data (type ''help <a href="matlab:help PlotCCG">PlotCCG</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CCG">CCG</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'full',
			full = lower(varargin{i+1});
			if ~isstring_FMAT(full,'on','off'),
				error('Incorrect value for property ''full'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
			full = strcmp(full,'on');
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).']);
	end
end

[~,nCols,nRows] = size(ccg); % reference series in columns, referenced series in rows
if nCols ~= nRows && ~full,
	warning('Unequal numbers of rows and columns - forcing ''full'' mode');
	full = true;
end

if full,
	for col = 1:nCols,
		for row = 1:nRows,
			subplot(nRows,nCols,col+nCols*(nRows-row));
			b = bar(t,ccg(:,col,row));
			xlim([min(t) max(t)]);
			set(b,'EdgeColor','k','FaceColor','k');
			if col == 1,
				ylabel(int2str(row));
			end
		end
		if col == row,
			xlabel('Time (s)');
		else
			set(gca,'xtick',[]);
		end
		if row == nRows,
			title(int2str(col));
		end
	end
else
	n = size(ccg,2);
	for col = 1:n,
		for row = col:n,
			subplot(n,n,col+n*(n-row));
			b = bar(t,ccg(:,col,row));
			xlim([min(t) max(t)]);
			if col == row,
				set(b,'EdgeColor','none','FaceColor','b');
			else
				set(b,'EdgeColor','k','FaceColor','k');
			end
			if col == 1,
				ylabel(int2str(row));
			end
		end
		if col == row,
			xlabel('Time (s)');
		else
			set(gca,'xtick',[]);
		end
		if row == n,
			title(int2str(col));
		end
	end
end


