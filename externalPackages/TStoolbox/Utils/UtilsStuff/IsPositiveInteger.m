%IsPositiveInteger - Test if array contains only positive integers.
%
%  USAGE
%
%    test = IsPositiveInteger(x,strict)
%
%    x              array to test
%    strict         optionally, test for strictly positive integers
%                   (default = 'off')


% Copyright (C) 2004-2006 by Michael Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

function test = IsPositiveInteger(x,strict)

if nargin == 1,
	strict = 'off';
end

x = x(:);

if strcmp(strict,'on'),
	test = isa(x,'numeric') & all(x==round(x)) & all(x>0);
else
	test = isa(x,'numeric') & all(x==round(x)) & all(x>=0);
end