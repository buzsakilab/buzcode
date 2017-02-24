function [data,eid,name] = DBGetValues(query,direction)

%DBGetValues - Group values for all variables that match given criteria.
%
%  Concatenate individual values for all variables matching the optional query.
%  Values must be scalars, vectors or matrices. Concatenation can be performed
%  horizontally or vertically. Matrices can also be concatenated along the
%  third dimension (z). By default, values of size MxN are concatenated
%  vertically if M<=N, or horizontally otherwise.
%
%  USAGE
%
%    [data,eid,name] = DBGetValues(query,direction)
%
%    query          optional list query (WHERE clause; see Example)
%    direction      optional concatenation direction ('h', 'v' or 'z')
%
%  OUTPUT
%
%    data           concatenation of individual values
%    eid            experiment ID (identifier string)
%    name           descriptive name (identifier string)
%
%  EXAMPLE
%
%    In this example, place cells were recorded on successive days while the
%    animal explored a maze, then during sleep. Firing fields were stored in a
%    database using eids like '20120213-Maze-(1,2)' (where 1,2 corresponds to
%    tetrode 1, cluster 2) and named 'FiringFields'. Mean firing rates during
%    sleep were named 'MeanRate'.
%
%    Get all firing fields (for details on query syntax, see an SQL manual):
%
%    [fields,eid1] = DBGetValues('eid like "%Maze%" and name="FiringField"');
%
%    Get all firing rates during sleep:
%
%    [rates,eid2] = DBGetValues('eid like "%Sleep%" and name="MeanRate"');
%
%  SEE
%
%    See also DBMatchValues, DBGetVariables, DBAddVariable, DBGetFigures.
%

% Copyright (C) 2007-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

data = [];
eid = {};
name = {};

% Parameters
if nargin < 1,
	query = '';
	direction = [];
elseif nargin == 1,
	if strcmp(lower(query),'h'),
		query = '';
		direction = 'h';
	elseif strcmp(lower(query),'v'),
		query = '';
		direction = 'v';
	elseif strcmp(lower(query),'z'),
		query = '';
		direction = 'z';
	else
		direction = [];
	end
end

% Query database
results = DBGetVariables(query,'output','full');
if isempty(results.v), return; end
eid = results.eid;
name = results.name;

first = results.v{1};
if ~isnumeric(first),
	error('Values are not numerical and cannot be concatenated (type ''help <a href="matlab:help DBGetValues">DBGetValues</a>'' for details).');
end

% Automatic concatenation direction
if isempty(direction),
	if size(first,1) <= size(first,2),
		direction = 'v';
	else
		direction = 'h';
	end
end

% Concatenation

% Determine sizes of retrieved values (arrays)
sizes = cellfun(@size,results.v,'UniformOutput',false);
% Make sure they all have the same number of dimensions
nDims = cellfun(@length,sizes);
if any(nDims~=nDims(1)),
	error('All values should have the same number of dimensions (type ''help <a href="matlab:help DBGetValues">DBGetValues</a>'' for details).');
end
% Determine 'circumscribed' array size
s = max(unique(vertcat(sizes{:}),'rows'),[],1);
% Extend all arrays
results.v = cellfun(@(x) ExtendArray(x,s),results.v,'uniformoutput',false);
% Concatenate in the appropriate direction
if strcmp(direction,'v'),
	data = vertcat(results.v{:});
elseif strcmp(direction,'h'),
	data = [results.v{:}];
else
	data = cat(3,results.v{:});
end
