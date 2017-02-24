function pairs = MatchPairs(list1,list2,varargin)

%MatchPairs - Pair nearest values in two lists.
%
%  USAGE
%
%    pairs = MatchPairs(list1,list2,<options>)
%
%    list1          first list
%    list2          second list
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'error'       maximum error (default = 0)
%    =========================================================================
%
%  OUTPUT
%
%    pairs        list of matched pairs (NaN indicates no match)
%
%  NOTE
%
%    This function uses list1 as the reference and is generally not symmetrical.






% ***************** TODO: Change PlotRippleStats to use Match instead of MatchPairs, then remove MatchPairs





% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
err = 0;

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help MatchPairs">MatchPairs</a>'' for details).');
end

% Check parameter sizes
if ~isdvector(list1),
	error('Parameter ''list1'' is not a vector (type ''help <a href="matlab:help MatchPairs">MatchPairs</a>'' for details).');
end
if ~isdvector(list2),
	error('Parameter ''list2'' is not a vector (type ''help <a href="matlab:help MatchPairs">MatchPairs</a>'' for details).');
end
if size(list2,1) == 1,
	list2 = list2';
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help MatchPairs">MatchPairs</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'error',
      err = varargin{i+1};
      if ~isdscalar(err,'>0'),
        error('Incorrect value for property ''error'' (type ''help <a href="matlab:help MatchPairs">MatchPairs</a>'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help MatchPairs">MatchPairs</a>'' for details).']);
  end
end

i = 1;
i2 = 1;
for i1 = 1:length(list1),
	% disp(['(' int2str(i1) ') ' num2str(list1(i1)) ' - (' int2str(i2) ') ' num2str(list2(i2)) ' ?']);
	while list1(i1) > list2(i2) + err,
		% disp([' (' int2str(i1) ') ' num2str(list1(i1)) ' > (' int2str(i2) ') ' num2str(list2(i2)) ' + ' num2str(err)]);
		pairs(i,:) = [NaN list2(i2)];
		i = i+1;
		i2 = i2+1;
	end
	if list1(i1) < list2(i2) - err,
		% disp([' (' int2str(i1) ') ' num2str(list1(i1)) ' < (' int2str(i2) ') ' num2str(list2(i2)) ' - ' num2str(err)]);
		pairs(i,:) = [list1(i1) NaN];
		i = i+1;
	else
		% disp([' (' int2str(i1) ') ' num2str(list1(i1)) ' ~ (' int2str(i2) ') ' num2str(list2(i2)) ' +- ' num2str(err)]);
		pairs(i,:) = [list1(i1) list2(i2)];
		i2 = i2+1;
		i = i+1;
	end
end
n = length(list2)-i2+1;
if n > 0,
	pairs = [pairs ; repmat(NaN,n,1) list2(i2:end)];
end