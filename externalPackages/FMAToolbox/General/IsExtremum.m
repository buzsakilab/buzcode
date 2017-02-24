function l = IsExtremum(samples,varargin)

%IsExtremum - Identify local maxima or minima.
%
% Find the maxima (resp. minima) in a signal by smoothing the signal and then
% finding the points where the derivative goes from positive to negative (resp.
% from negative to positive).
%
%  USAGE
%
%    status = IsExtremum(samples,<options>)
%
%    samples        an Nx(M+1) matrix of (timestamp,value1,...,valueM) t-uples
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'mode'        'maxima' or 'minima' (default = 'maxima')
%     'type'        'linear' if samples are linear values (default),
%                   'circular' otherwise
%     'smooth'      standard deviation for Gaussian kernel (0 = no smoothing,
%                   default)
%    =========================================================================
%
%  OUTPUT
%
%    status         an NxM logical matrix, i.e. status(i,j)==1 if samples(i,j)
%                   is an extremum
%

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
mode = 'maxima';
smooth = 0;
type = 'linear';

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help IsExtremum">IsExtremum</a>'' for details).');
end

% Parse parameter list
for j = 1:2:length(varargin),
	if ~ischar(varargin{j}),
		error(['Parameter ' num2str(j+7) ' is not a property (type ''help <a href="matlab:help IsExtremum">IsExtremum</a>'' for details).']);
	end
	switch(lower(varargin{j})),
		case 'mode',
			mode = lower(varargin{j+1});
			if ~isstring_FMAT(mode,'minima','maxima'),
				error('Incorrect value for ''mode'' (type ''help <a href="matlab:help IsExtremum">IsExtremum</a>'' for details).');
			end
		case 'type',
			type = lower(varargin{j+1});
			if ~isstring_FMAT(type,'linear','circular'),
				error('Incorrect value for ''type'' (type ''help <a href="matlab:help Diff">Diff</a>'' for details).');
			end
		case 'smooth',
			smooth = varargin{j+1};
			if ~isdvector(smooth,'>=0') || length(smooth) > 2,
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help IsExtremum">IsExtremum</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{j}) ''' (type ''help <a href="matlab:help IsExtremum">IsExtremum</a>'' for details).']);
	end
end

if size(samples,2) ~= 2,
	smooth = [smooth 0];
end

derivative = Diff(samples,'smooth',smooth,'type',type);

previous = derivative(1:end-1,2:end);
next = derivative(2:end,2:end);

if strcmp(mode,'maxima'),
	l = previous > 0 & next < 0;
else
	l = previous < 0 & next > 0;
end
l = [l;zeros(1,size(l,2))];
l = logical(l);

