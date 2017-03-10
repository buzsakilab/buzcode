function [spectrum,f,s] = MTPointSpectrum(times,varargin)

%MTPointSpectrum - Compute point process spectrum by multi-taper estimation.
%
%  The spectrum is computed using the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  USAGE
%
%    [spectrum,f,s] = MTPointSpectrum(times,<options>)
%
%    times          event times, e.g. spike times
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   sampling rate (in Hz, default = 20,000 Hz)
%     'range'       frequency range (in Hz) (default = all)
%     'tapers'      relative resolution and order of the tapers [NW K]
%                   (default = [3 5])
%     'pad'         FFT padding (see help for <a href="matlab:help mtspectrumpt">mtspectrumpt</a>) (default = 0)
%     'show'        plot spectrum (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    spectrum       spectrum
%    s              error
%    f              frequency bins
%
%  DEPENDENCIES
%
%    This function requires the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  SEE
%
%    See also MTSpectrum, MTSpectrogram, SpectrogramBands.

% Copyright (C) 2010-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure chronux is installed and functional
CheckChronux('mtspectrumpt');

% Defaults
frequency = 20000;
range = [];
show = 'off';
tapers = [3 5];
pad = 0;
err = [1 0.95];

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help MTPointSpectrum">MTPointSpectrum</a>'' for details).');
end

% Check parameter sizes
if size(times,2) ~= 1 && size(times,2) ~= 2,
	error('Parameter ''times'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTPointSpectrum">MTPointSpectrum</a>'' for details).');
end

% Parse parameter list
v = {};
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help MTPointSpectrum">MTPointSpectrum</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'frequency',
			frequency = varargin{i+1};
			if ~isdscalar(frequency,'>0'),
				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help MTPointSpectrum">MTPointSpectrum</a>'' for details).');
			end
		case 'range',
			range = varargin{i+1};
			if ~isdvector(range,'#2','<','>=0'),
				error('Incorrect value for property ''range'' (type ''help <a href="matlab:help MTPointSpectrum">MTPointSpectrum</a>'' for details).');
			end
		case 'tapers',
			tapers = varargin{i+1};
			if ~isivector(tapers,'#2'),
				error('Incorrect value for property ''tapers'' (type ''help <a href="matlab:help MTPointSpectrum">MTPointSpectrum</a>'' for details).');
			end
		case 'pad',
			pad = varargin{i+1};
			if ~isiscalar(pad,'>-1'),
				error('Incorrect value for property ''pad'' (type ''help <a href="matlab:help MTPointSpectrum">MTPointSpectrum</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring_FMAT(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help MTPointSpectrum">MTPointSpectrum</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help MTPointSpectrum">MTPointSpectrum</a>'' for details).']);
	end
	if ~strcmp(varargin{i},'show'), v = {v{:},varargin{i:i+1}}; end
end

% Compute point spectrogram and moments
[spectrogram,~,f] = MTPointSpectrogram(lfp,v{:});
spectrogram = spectrogram';
spectrum = mean(spectrogram);
s = var(spectrogram);

% Plot log, i.e. mean E[log(spectrogram)] and stdev sqrt(Var[log(spectrogram)])
% (see http://en.wikipedia.org/wiki/Taylor_expansions_for_the_moments_of_functions_of_random_variables)
if strcmp(lower(show),'on'),
	figure;
	%logSpectrum = log(spectrum);         % classic value
	logSpectrum = log(mu)-v./(2*mu.*mu);  % corrected value, valid for mu/s > 1.5
	logS = sqrt(v./mu.^2);               % valid for mu/s > 2.5
	PlotMean(f,logSpectrum,logSpectrum-logS,logSpectrum+logS,':');
	xlabel('Frequency (Hz)');
	ylabel('Power');
	title('Power Spectrum');
end
