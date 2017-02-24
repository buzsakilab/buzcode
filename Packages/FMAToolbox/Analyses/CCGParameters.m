function [times,id,group] = CCGParameters(varargin)

%CCGParameters - Reformat time series for CCG computation
%
%  USAGE
%
%    [times,id,group] = CCGParameters(series1,group1,series2,group2,...)
%
%    series1...     time series (one column) with optional IDs (in the second
%                   column, such as obtained from <a href="matlab:help GetSpikes">GetSpikes</a>, see below)
%    group1...      optional group number for each series
%
%  SEE
%
%    See also CCG.

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CCGParameters">CCGParameters</a>'' for details).');
end

data = [];

% Parse parameter list
i = 1;
while i <= length(varargin),

	% Time series i (and optional IDs)
	times = varargin{i};
	i = i + 1;
	if ~isdvector(times) && ~(isdmatrix(times) && size(times,2) == 2),
		error(['Incorrect time series at parameter #' int2str(i) ' (type ''help <a href="matlab:help CCGParameters">CCGParameters</a>'' for details).']);
	end
	if size(times,2) > 1,
		id = times(:,2);
		if ~isivector(id,'>0'),
			error(['Incorrect time series at parameter #' int2str(i) ' (type ''help <a href="matlab:help CCGParameters">CCGParameters</a>'' for details).']);
		end
	else
		id = ones(size(times(:,1)));
	end
	if length(unique(id)) ~= max(id),
		error('Incorrect IDs (type ''help <a href="matlab:help CCGParameters">CCGParameters</a>'' for details).');
	end
	times = times(:,1);
	
	% Groups i
	if i > length(varargin),
		group = ones(size(id));
	else
		group = varargin{i};
		% Is this a group, or the next time series?
		if ~isivector(group) || (any(group>1) && any(group>2)),
			group = ones(size(id));
		else
			if length(group) == 1,
				group = repmat(group,size(id),1);
			elseif length(group) ~= length(id),
				error(['Incorrect groups at parameter #' int2str(i) ' (type ''help <a href="matlab:help CCGParameters">CCGParameters</a>'' for details).']);
			end
			i = i + 1;
		end
	end
	group = group(:);
	
	% Concatenate with previous times series, IDs and groups
	if ~isempty(data), id = id+max(data(:,2)); end
	data = [data ; times id group];

end

times = data(:,1);
id = data(:,2);
group = data(:,3);
