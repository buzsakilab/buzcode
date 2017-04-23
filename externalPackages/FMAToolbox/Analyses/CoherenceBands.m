function bands = CoherenceBands(coherogram,frequencies,varargin)

%CoherenceBands - Determine running coherence in physiological bands.
%
%  USAGE
%
%    bands = CoherenceBands(coherogram,frequencies,<options>)
%
%    coherogram     coherogram obtained using <a href="matlab:help MTCoherogram">MTCoherogram</a>
%    frequencies    frequency bins for coherogram
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      smoothing for ratio (0 = no smoothing) (default = 2)
%     'custom'      set custom band (default= [0 250])
%     'delta'       set delta band (default = [0 4])
%     'theta'       set theta band (default = [7 10])
%     'spindles'    set spindle band (default = [10 20])
%     'lowGamma'    set low gamma band (default = [30 80])
%     'highGamma'   set high gamma band (default = [80 120])
%     'ripples'     set ripple band (default = [100 250])
%    =========================================================================
%
%  OUTPUT
%
%    bands.custom     coherence in the custom band
%    bands.delta      coherence in the delta band
%    bands.theta      coherence in the theta band
%    bands.spindles   coherence in the spindle band
%    bands.lowGamma   coherence in the low gamma band
%    bands.highGamma  coherence in the high gamma band
%    bands.ripples    coherence in the ripple band
%
%  SEE
%
%    See also MTCoherogram.

% Copyright (C) 2013 by Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
smooth = 2;
custom = [0 250];
delta = [0 4];
theta = [7 10];
spindles = [10 20];
lowGamma = [30 80];
highGamma = [80 120];
ripples = [100 250];

% Check number of parameters
if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).');
end

% Check parameter sizes
if ~isdmatrix(coherogram),
	error('Parameter ''coherogram'' is not a matrix (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).');
end
if ~isdvector(frequencies) | length(frequencies) ~= size(coherogram,1),
	error('Parameter ''frequencies'' is not a vector or does not match coherogram size (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'smooth',
			smooth = varargin{i+1};
			if ~isdscalar(smooth,'>=0'),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).');
			end
		case 'custom',
			custom = varargin{i+1};
			if ~isdvector(delta,'>=0','<','#2'),
				error('Incorrect value for property ''delta'' (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).');
			end
		case 'delta',
			delta = varargin{i+1};
			if ~isdvector(delta,'>=0','<','#2'),
				error('Incorrect value for property ''delta'' (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).');
			end
		case 'theta',
			theta = varargin{i+1};
			if ~isdvector(theta,'>=0','<','#2'),
				error('Incorrect value for property ''theta'' (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).');
			end
		case 'spindles',
			spindles = varargin{i+1};
			if ~isdvector(spindles,'>=0','<','#2'),
				error('Incorrect value for property ''spindles'' (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).');
			end
		case 'lowgamma',
			lowGamma = varargin{i+1};
			if ~isdvector(lowGamma,'>=0','<','#2'),
				error('Incorrect value for property ''lowGamma'' (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).');
			end
		case 'highgamma',
			highGamma = varargin{i+1};
			if ~isdvector(highGamma,'>=0','<','#2'),
				error('Incorrect value for property ''highGamma'' (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).');
			end
		case 'ripples',
			ripples = varargin{i+1};
			if ~isdvector(ripples,'>=0','<','#2'),
				error('Incorrect value for property ''ripples'' (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SpectrogramBands">SpectrogramBands</a>'' for details).']);
	end
end

% Select relevant frequency bins
customBins = frequencies >= custom(1) & frequencies <= custom(2);
thetaBins = frequencies >= theta(1) & frequencies <= theta(2);
deltaBins = frequencies >= delta(1) & frequencies <= delta(2);
spindleBins = frequencies >= spindles(1) & frequencies <= spindles(2);
lowGammaBins = frequencies >= lowGamma(1) & frequencies <= lowGamma(2);
highGammaBins = frequencies >= highGamma(1) & frequencies <= highGamma(2);
rippleBins = frequencies >= ripples(1) & frequencies <= ripples(2);

% Compute physiological bands
bands.custom = mean(coherogram(customBins,:))';
bands.theta = mean(coherogram(thetaBins,:))';
bands.delta = mean(coherogram(deltaBins,:))';
bands.spindles = mean(coherogram(spindleBins,:))';
bands.lowGamma = mean(coherogram(lowGammaBins,:))';
bands.highGamma = mean(coherogram(highGammaBins,:))';
bands.ripples = mean(coherogram(rippleBins,:))';

