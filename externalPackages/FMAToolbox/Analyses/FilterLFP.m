function filtered = FilterLFP(lfp,varargin)

%FilterLFP - Filter the local field potentials, e.g. in the theta band.
%
% Filter the local field potentials using a cheby2 filter.
%
%  USAGE
%
%    filtered = FilterLFP(lfp,<options>)
%
%    lfp            local field potentials, list of (t,v1,v2,v3...) tuples
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'passband'    frequency range (default = [4 10], 'theta')
%     'order'       filter order (default = 4)
%     'ripple'      filter ripple (default = 20)
%     'nyquist'     nyquist frequency (default = 625)
%     'filter'      choose filter type between 'cheby2' (default) and 'fir1'
%    =========================================================================
%
%  NOTE
%
%    The passband can be supplied either explicitly, e.g. [30 80], or by name,
%    by choosing among the following predefined frequency bands:
%
%        delta      [0 4]
%        theta      [4 10]
%        spindles   [10 20]
%        gamma      [30 80]
%        ripples    [100 250]

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

error('this function is now deprecated, try using bz_Filter.')

% Default values
passband = [4 10];

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help FilterLFP">FilterLFP</a>'' for details).');
end

% Check parameter sizes
if ~isdmatrix(lfp),
	error('Parameter ''lfp'' is not a matrix (type ''help <a href="matlab:help FilterLFP">FilterLFP</a>'' for details).');
end

% Parse parameter list
i = 1;
while i <= length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help FilterLFP">FilterLFP</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'passband',
			passband = varargin{i+1};
			if isstring_FMAT(passband),
				switch(lower(passband)),
					case 'delta',
						passband = [0 4];
					case 'theta',
						passband = [4 10];
					case 'spindles',
						passband = [10 20];
					case 'gamma',
						passband = [30 80];
					case 'ripples',
						passband = [100 250];
					otherwise,
						error(['Unknown frequency band ''' passband ''' (type ''help <a href="matlab:help FilterLFP">FilterLFP</a>'' for details).']);
				end
			elseif ~isdvector(passband,'#2','>=0'),
				error('Incorrect value for ''passband'' (type ''help <a href="matlab:help FilterLFP">FilterLFP</a>'' for details).');
			end
			varargin = {varargin{[1:(i-1) (i+2):length(varargin)]}};
		otherwise,
			i = i+2;
	end
end
filtered = Filter(lfp,'passband',passband,varargin{:});