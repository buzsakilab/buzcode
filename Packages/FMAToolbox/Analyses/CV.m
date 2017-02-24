function [coeff,local] = CV(timestamps,varargin)

%CV - Compute coefficient of variation for a point process.
%
% Three alternative measures can be computed: regular CV, CV in operational time
% (Gestri and Petracchi, 1970; Nawrot et al., 2008) or CV2 (Holt et al., 2006).
%
% While the CV is typically computed for consecutive events, it can also be
% applied to intervals between events i and i+2 (order 2), or between events
% i and i+3 (order 3), etc. to study higher order variability.
%
% Although operational time is a useful concept, in practice the estimate of CV
% strongly depends on the estimated instantaneous frequency of the point process
% (here, instantaneous frequency is estimated using an adaptive kernel filtering
% method, see <a href="matlab:help Frequency">Frequency</a>). Use this measure with care.
%
%  USAGE
%
%    [coeff,local] = CV(timestamps,<options>)
%
%    timestamps     point process
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'measure'     either 'cv' (default), 'cvo' (operational time) or 'cv2'
%     'order'       interval order (default = 1)
%     'method'      to estimate instantaneous frequency (default = 'fixed')
%     'binSize'     to estimate instantaneous frequency (default = 0.001 s)
%     'smooth'      to estimate instantaneous frequency (default = 25)
%    =========================================================================
%
%  OUTPUT
%
%    coeff          (mean) coefficient of variation
%    local          time-dependent coefficient of variation (CV2 only)
%

% Copyright (C) 2010-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
measure = 'cv';
binSize = 0.001;
smooth = 25;
order = 1;
method = 'fixed';

% Check parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CV">CV</a>'' for details).');
end
if ~isdvector(timestamps),
  error('Incorrect timestamps (type ''help <a href="matlab:help CV">CV</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help CV">CV</a>'' for details).']);
  end
  switch(lower(varargin{i})),
		case 'measure',
			measure = lower(varargin{i+1});
			if ~isstring_FMAT(measure,'cv','cv2','cvo'),
				error('Incorrect value for property ''measure'' (type ''help <a href="matlab:help CV">CV</a>'' for details).');
			end
		case 'method',
			method = lower(varargin{i+1});
			if ~isstring_FMAT(method,'fixed','adaptive','inverse'),
				error('Incorrect value for property ''method'' (type ''help <a href="matlab:help CV">CV</a>'' for details).');
			end
		case 'order',
			order = varargin{i+1};
			if ~isiscalar(order,'>0'),
				error('Incorrect value for property ''order'' (type ''help <a href="matlab:help CV">CV</a>'' for details).');
			end
		case 'binsize',
			binSize = varargin{i+1};
			if ~isdscalar(binSize,'>0'),
				error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help CV">CV</a>'' for details).');
			end
		case 'smooth',
			smooth = varargin{i+1};
			if ~isdscalar(smooth,'>=0'),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help CV">CV</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help CV">CV</a>'' for details).']);
  end
end

switch(measure),
	case 'cv',
		dt = ndiff(timestamps,order);
		coeff = std(dt)/mean(dt);
		local = [];
	case 'cvo',
		% Operational time
		% 1) Compute instantaneous frequency by adaptive filtering
		frequency = Frequency(timestamps,'binSize',binSize,'smooth',smooth,'method',method);
		% 2) Integrate frequency to yield operational time
		operational = [frequency(:,1) cumsum(frequency(:,2))*binSize];
		operational = Interpolate(operational,timestamps);
		operational = operational(:,2)+timestamps(1)-operational(1,2);
		dto = ndiff(operational,order);
		coeff = std(dto)/mean(dto);
		local = [];
	case 'cv2',
		dt = ndiff(timestamps,order);
		dt1 = dt(1:end-1);
		dt2 = dt(2:end);
		x = dt1./dt2; % ratio of consecutive iter-event intervals
		local = 2*abs(x-1)./(x+1);
		coeff = mean(local);
end

function y = ndiff(x,n)

if n == 1,
	y = diff(x);
else
	y = x((n+1):end)-x(1:end-n);
end