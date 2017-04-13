%GetNextField - Get next field from a Batch object (iterator mode).
%
% This is a helper class to easily read and parse batch files.
%
%  USAGE
%
%    [b,field] = GetNextField(b);
%
%    b              batch object

% Copyright (C) 2007 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [b,field] = GetNextField(b)

% Check number of parameters
if nargin ~= 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help GetNextField">GetNextField</a>'' for details).');
end

b.currentField = b.currentField+1;
if b.currentItem < 1 || b.currentItem > size(b.field,1) || b.currentField > size(b.field,2),
	field = [];
else
	field = b.field{b.currentItem,b.currentField};
end
