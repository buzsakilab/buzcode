%minmax - Show min and max values for any number of variables.
%
%  USAGE
%
%    minmax(x,y,...)
%
%    x,y,...        variables
%

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function minmax(varargin)

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help minmax">minmax</a>'' for details).');
end

for i = 1:length(varargin),
	str = '';
	% Variable index (min)
	str = [str '   [' int2str(i) '-] '];
	% Variable values (min)
	s = min(varargin{i});
	for j = 1:length(s),
		str = [str num2str(s(j)) ' '];
	end
	str(end) = [];
	disp(str);
	str = '';
	% Variable index (max)
	str = [str '   [' int2str(i) '+] '];
	% Variable values (max)
	s = max(varargin{i});
	for j = 1:length(s),
		str = [str num2str(s(j)) ' '];
	end
	str(end) = [];
	disp(str);
end
