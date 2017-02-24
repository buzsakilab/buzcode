%sz - Show sizes for any number of variables.
%
%  USAGE
%
%    sz(x,y,...)
%
%    x,y,...        variables
%

% Copyright (C) 2011-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function sz(varargin)

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help sz">sz</a>'' for details).');
end

for i = 1:length(varargin),
	str = '';
	s = size(varargin{i});
	% Variable index
	str = [str '   [' int2str(i) '] '];
	% Variable size
	for j = 1:length(s),
		str = [str int2str(s(j)) 'x'];
	end
	str(end) = [];
	disp(str);
end
