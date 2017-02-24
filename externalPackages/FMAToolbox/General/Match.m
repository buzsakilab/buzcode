function [matched,only1,only2] = Match(list1,list2,varargin)

%Match - Replace values in one list with closest values in a second list.
%
%  USAGE
%
%    [matched,only1,only2] = Match(list1,list2,<options>)
%
%    list1          first list
%    list2          second list
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'match'       optional matching criterion: 'up', 'down' (default) or
%                   'closest' (see below)
%     'error'       matched values for which the difference exceeds this
%                   threshold are set to NaN (default = inf)
%    =========================================================================
%
%  OUTPUT
%
%    matched        values in list2 that matched those of list1
%    only1          elements in list1 that were not matched by any element in list2
%    only2          elements in list2 that did not match any element in list1
%
%    only1 and only2 are lists of logical indices.
%
%  NOTE
%
%    Pairs are matched in the following way. Assume the lists are sorted in
%    ascending order and let list1 = [a(1)...a(n)] and list2 = [b(1)...b(m)].
%    We are trying to match each a(i) with one of the b(1..m). Suppose
%
%       b(1) <= ... <= b(j) <= a(i) <= b(j+1) <= ... <= b(m)
%
%    Then a(i) will be matched with either b(j) or b(j+1), depending on the
%    'match' option:
%
%       1) 'down' (default) => b(j)
%       2) 'up' => b(j+1)
%       3) 'closest' => the one that minimizes the distance
%          (b(j) if the distances are equal)
%
%    If however the distance between a(i) and b(j) (or b(j+1)) is too large
%    (see option 'error'), then none of the b(1..m) matches a(i).
%
%    Note that this function is generally asymmetrical.

% Copyright (C) 2004-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
match = 'down';
err = inf;

% Check parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Match">Match</a>'' for details).');
end
if ~isdvector(list1) | ~isdvector(list2),
	error('Incorrect sizes: ''list1'' and ''list2'' must be vectors (type ''help <a href="matlab:help Match">Match</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Match">Match</a>'' for details).']);
  end
  switch(lower(varargin{i})),

    case 'match',
		match = lower(varargin{i+1});
		if ~isstring_FMAT(match,'up','down','closest'),
			error('Incorrect value for property ''match'' (type ''help <a href="matlab:help Match">Match</a>'' for details).');
      end

	case 'error',
		err = varargin{i+1};
		if ~isdscalar(err,'>=0'),
			error('Incorrect value for property ''error'' (type ''help <a href="matlab:help Match">Match</a>'' for details).');
		end

    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Match">Match</a>'' for details).']);

  end
end

% Sort elements in both list so that list1 = [a1 a2 a3 ... an] and list2 = [b1 b2 b3 ... bm]
[list1,idx1] = sort(list1(:));
[list2,idx2] = sort(list2(:));

% Find 'up' and 'down' matching indices
% (We need both, because this will allow us to compare matches if the criterion is 'closest')
n = length(list2);
up = MatchUpIndices(list1,list2);
down = up-1;

% List corresponding values in list2
% Special cases: when matching down (resp. up), values in list1 lesser (resp. greater)
% than all values in list2 cannot be matched; set them to NaN
goodUp = up>0 & up<=n;
matchedUp = nan(size(up));
matchedUp(goodUp) = list2(up(goodUp));
goodDown = down>0 & down<=n;
matchedDown = nan(size(down));
matchedDown(goodDown) = list2(down(goodDown));

% Use match criterion
switch(match),
	case 'up',
		matched = matchedUp;
	case 'down',
		matched = matchedDown;
	case 'closest',
		[unused,idx] = min(abs([matchedUp matchedDown]-[list1 list1]),[],2);
		matched = matchedUp;
		matched(idx==2) = matchedDown(idx==2);
end

% Check error threshold
bad = abs(matched-list1) > err+eps; % eps is to compensate for inevitable numerical rounding errors such as 1.1-1.0>0.1
matched(bad) = nan;

% Output
only1 = isnan(matched);
matched = matched(~only1);
only1(idx1) = only1;
only2 = logical(zeros(size(list2)));
m = unique(matched);
[unused,i] = setdiff(list2,m);
only2(i) = 1;
only2(idx2) = only2;

