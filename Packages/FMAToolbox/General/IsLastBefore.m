function l = IsLastBefore(timestamps,list,varargin)

%IsLastBefore - Identify last item before each of a list of timestamps.
%
% For each element i of a list of test timestamps, find in a list of reference timestamps
% the last element j smaller than i. Returns a list of logical indices (the elements of the
% test list which are smaller than all the elements of the reference timestamps are ignored).
%
%  USAGE
%
%    logical = IsLastBefore(timestamps,list,<options>)
%
%    timestamps     a list of reference timestamps
%    list           a list of test timestamps
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'strict'      'on' for strict (<) comparisons, 'off' otherwise (<=)
%                   (default = 'off')
%    =========================================================================
%
%  SEE
%
%    See also IsFirstAfter.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

strict = 'off';

if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help IsLastBefore">IsLastBefore</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help IsLastBefore">IsLastBefore</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'strict',
		strict = lower(varargin{i+1});
		if ~isstring_FMAT(strict,'on','off'),
			error('Incorrect value for property ''strict'' (type ''help <a href="matlab:help IsLastBefore">IsLastBefore</a>'' for details).');
		end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help IsLastBefore">IsLastBefore</a>'' for details).']);
  end
end

% Make sure 'timestamps' looks like [t1;t2;...;tn] and 'list' like [l1 l2 ... lm]
if min(size(timestamps)) ~= 1 | min(size(list)) ~= 1,
	error('Incorrect sizes (must be vectors)');
end
if size(timestamps,1) == 1,
	timestamps = timestamps';
end
if size(list,2) == 1,
	list = list';
end

% Now compare [t1 t1 ... t1;t2 t2 ... t2;...;tn tn ... tn] with [l1 l2 ... lm;l1 l2 ... lm;...;l1 l2 ... lm]
if strcmp(strict,'on'),
	where = repmat(timestamps,1,length(list)) < repmat(list,length(timestamps),1);
else
	where = repmat(timestamps,1,length(list)) <= repmat(list,length(timestamps),1);
end
% This yielded something like [1 1 1;1 1 1;...;0 1 1;...;0 0 1;...;0 0 0;...;0 0 0]
% where the transition from 1 to 0 in each column indicates the index of the last element
% in 'timestamps' inferior to the corresponding element in 'list'. Now, find these transitions.
% (we append ones before and zeros after 'where' to handle the special cases where the first
% element in 'timestamps' inferior to the element of 'list' is the last element of 'timestamps'
% or does not exist).
where = diff([ones(1,size(where,2));where;zeros(1,size(where,2))]);
% Convert the resulting logical matrix into a list of logical indices
l = sum(where,2) == -1;
l = l(2:length(timestamps)+1);
l = logical(l);
