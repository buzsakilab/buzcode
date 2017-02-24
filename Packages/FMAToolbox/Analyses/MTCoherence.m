function [coherence,f,s] = MTCoherence (lfp1,lfp2,varargin)

%MTCoherence - Compute Coherence by multi-taper estimation.
%
%  The spectrum is computed as the average spectrogram, using the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  USAGE
%
%    [coherence,f,s] = MTSpectrum(lfp1,lfp2,<options>)
%
%    lfp1            wide-band LFP 1 (one channel).
%    lfp2            wide-band LFP 2 (one channel).
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   sampling rate (in Hz) (default = from timestamps if
%                   available, otherwise 1250Hz)
%     'range'       frequency range (in Hz) (default = all)
%     'window'      duration (in s) of the time window (default = 5)
%     'overlap'     overlap between successive windows (default = window/2)
%     'step'        step between successive windows (default = window/2)
%     'tapers'      relative resolution and order of the tapers [NW K]
%                   (default = [3 5])
%     'pad'         FFT padding (see help for <a href="matlab:help mtspecgramc">mtspectrumc</a>) (default = 0)
%     'error'       error type (see help for <a href="matlab:help mtspecgramc">mtspectrumc</a>) (default = [1 0.95])
%     'show'        plot log spectrum (default = 'off')
%    =========================================================================
%
%  NOTES
%
%    The LFP can be provided either as a time stamped matrix (list of time-voltage
%    pairs), or as a voltage vector - in which case the frequency must be specified.
%
%    The time displacement between successive short time spectra can be supplied
%    either as a 'step' (explicit time difference) or as an 'overlap' (between
%    successive time windows).
%
%  OUTPUT
%
%    coherence      coherence spectrum
%    s              error
%    f              frequency bins
%
%  DEPENDENCIES
%
%    This function requires the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  SEE
%
%    See also MTCoherogram, SpectrogramBands.

% Nov 2013, Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure chronux is installed and functional
CheckChronux;

% Defaults
f = 1250;
frequency = [];
window = 5;
range = [];
overlap = [];
step = [];
show = 'off';
tapers = [3 5];
pad = 0;
err = [1 0.95];

% Check dependencies
if isempty(which('mtspecgramc')),
	error('This function requires the <a href="http://www.chronux.org">chronux</a> toolbox by P. Mitra, which does not appear to be installed on this system.');
end

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
end

% Check parameter sizes
if size(lfp1,2) ~= 1 && size(lfp1,2) ~= 2,
	error('Parameter ''lfp'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
end
if size(lfp2,2) ~= 1 && size(lfp2,2) ~= 2,
	error('Parameter ''lfp'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
end


% Parse parameter list
v = {};
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'frequency',
			frequency = varargin{i+1};
			if ~isdscalar(frequency,'>0'),
				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
			end
		case 'range',
			range = varargin{i+1};
			if ~isdvector(range,'#2','<','>=0'),
				error('Incorrect value for property ''range'' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
			end
		case 'window',
			window = varargin{i+1};
			if ~isdscalar(window,'>0'),
				error('Incorrect value for property ''window'' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
			end
		case 'overlap',
			overlap = varargin{i+1};
			if ~isdscalar(overlap,'>0'),
				error('Incorrect value for property ''overlap'' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
			end
		case 'step',
			step = varargin{i+1};
			if ~isdscalar(step,'>0'),
				error('Incorrect value for property ''step'' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
			end
		case 'tapers',
			tapers = varargin{i+1};
			if ~isivector(tapers,'#2','>0'),
				error('Incorrect value for property ''tapers'' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
			end
		case 'pad',
			pad = varargin{i+1};
			if ~isiscalar(pad,'>-1'),
				error('Incorrect value for property ''pad'' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).']);
	end
	if ~strcmp(varargin{i},'show'), v = {v{:},varargin{i:i+1}}; end
end

% Compute coherogram and moments
[coherogram,t,f] = MTCoherogram(lfp1,lfp2,v{:});
coherogram = coherogram';
mu = mean(coherogram);
v = var(coherogram);

% Compute coherence
coherence = mu;
s = sqrt(v);

% Plot log, i.e. mean E[log(spectrogram)] and stdev sqrt(Var[log(spectrogram)])
% (see http://en.wikipedia.org/wiki/Taylor_expansions_for_the_moments_of_functions_of_random_variables)
if strcmp(lower(show),'on'),
	figure;
	PlotMean(f,coherence,coherence-s,coherence+s,':');
	ylim([0 1]);
end
