function m = Bright(n,varargin)

%Bright - Bright colormap (similar to HSV or JET, but brighter).
%
%  USAGE
%
%    m = Bright(n,<options>)
%
%    n              optional number of rows in output matrix (default = 100)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'gamma'       gamma correction (1 = no correction, default)
%     'type'        either 'linear' or 'circular' (default 'linear')
%    =========================================================================

% Copyright (C) 2009-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
gamma = 1;
type = 'linear';

% Optional parameter
if nargin < 1,
	n = 100;
elseif ischar(n),
	varargin = {n,varargin{:}};
	n = 100;
elseif ~isdscalar(n,'>=0'),
	error('Incorrect value for ''n'' (type ''help <a href="matlab:help Bright">Bright</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Bright">Bright</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'gamma',
      gamma = varargin{i+1};
      if ~isdscalar(gamma,'>=0'),
        error('Incorrect value for property ''gamma'' (type ''help <a href="matlab:help Bright">Bright</a>'' for details).');
      end
    case 'type',
    	type = lower(varargin{i+1});
      if ~isstring_FMAT(type,'linear','circular'),
        error('Incorrect value for property ''type'' (type ''help <a href="matlab:help Bright">Bright</a>'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Bright">Bright</a>'' for details).']);
  end
end

if strcmp(type,'linear'),
	hsv(:,1) = (2/3) - (2/3)*linspace(0,1,n)'.^gamma;
else
	hsv(:,1) = linspace(0,1,n)'.^gamma;
end
hsv(:,2) = 1;
hsv(:,3) = 1;
m = hsv2rgb(hsv);
