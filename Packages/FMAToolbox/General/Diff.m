function derivative = Diff(samples,varargin)

%Diff - Differentiate.
%
%  USAGE
%
%    derivative = Diff(samples,<options>)
%
%    samples        data to differentiate
%                   measured in number of samples (default = no smoothing)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'type'        'linear' if samples are linear values (default),
%                   'circular' otherwise
%     'smooth'      standard deviation for Gaussian kernel (default = 0)
%    =========================================================================
%
%  SEE
%
%    See also Interpolate.
%

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
type = 'linear';
smooth = 0;

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Diff">Diff</a>'' for details).');
end

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Diff">Diff</a>'' for details).');
end

% Parse parameter list
for j = 1:2:length(varargin),
	if ~ischar(varargin{j}),
		error(['Parameter ' num2str(j+7) ' is not a property (type ''help <a href="matlab:help Diff">Diff</a>'' for details).']);
	end
	switch(lower(varargin{j})),
		case 'type',
			type = lower(varargin{j+1});
			if ~isstring_FMAT(type,'linear','circular'),
				error('Incorrect value for ''type'' (type ''help <a href="matlab:help Diff">Diff</a>'' for details).');
			end
		case 'smooth',
			smooth = varargin{j+1};
			if ~isdvector(smooth,'>=0') || length(smooth) > 2,
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help Diff">Diff</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{j}) ''' (type ''help <a href="matlab:help Diff">Diff</a>'' for details).']);
	end
end

% Circular?
if strcmp(type,'circular'),
	samples(:,2:end) = unwrap(samples(:,2:end));
end

% Smooth
if smooth ~= 0,
	if length(samples(1,2:end)) >= 2,
		smooth = [smooth 0];
	end
	samples(:,2:end) = Smooth(samples(:,2:end),smooth);
end

% Differentiate
derivative = zeros(size(samples));
derivative(:,1) = samples(:,1);
t = samples(:,1);
dt = diff(t);
t_diff = t(1:end-1)+dt/2;
for i = 2:size(samples,2),
	d0 = diff(samples(:,i))./dt;
	d1 = interp1(t_diff,d0,t(2:end-1,1));
	derivative(:,i) = [d0(1);d1;d0(end)];
end
