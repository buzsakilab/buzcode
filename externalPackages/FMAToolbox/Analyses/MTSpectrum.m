function [spectrum,f,s] = MTSpectrum(lfp,varargin)

%MTSpectrum - Compute LFP spectrum by multi-taper estimation.
%
%  The spectrum is computed using the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  USAGE
%
%    [spectrum,f,s] = MTSpectrum(lfp,<options>)
%
%    lfp            wide-band LFP (one channel).
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   sampling rate (in Hz) (default = from timestamps if
%                   available, otherwise 1250Hz)
%     'range'       frequency range (in Hz) (default = all)
%     'tapers'      relative resolution and order of the tapers [NW K]
%                   (default = [3 5])
%     'pad'         FFT padding (see help for <a href="matlab:help mtspecgramc">mtspectrumc</a>) (default = 0)
%     'error'       error type (see help for <a href="matlab:help mtspecgramc">mtspectrumc</a>) (default = [1 0.95])
%     'show'        plot spectrum (default = 'off')
%    =========================================================================
%
%  NOTE
%
%    The lfp can be provided either as a time stamped matrix (list of time-voltage
%    pairs), or as a voltage vector - in which case the frequency must be specified.
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
%    See also MTSpectrogram, SpectrogramBands.

% Copyright (C) 2010-2011 by Michaël Zugaro
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
range = [];
show = 'off';
tapers = [3 5];
pad = 0;
err = [1 0.95];

% Check dependencies
if isempty(which('mtspectrumc')),
	error('This function requires the <a href="http://www.chronux.org">chronux</a> toolbox by P. Mitra, which does not appear to be installed on this system.');
end

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
end

% Check parameter sizes
if size(lfp,2) ~= 1 && size(lfp,2) ~= 2,
	error('Parameter ''lfp'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
end

% Parse parameter list
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
		case 'tapers',
			tapers = varargin{i+1};
			if ~isivector(tapers,'#2'),
				error('Incorrect value for property ''tapers'' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
			end
		case 'pad',
			pad = varargin{i+1};
			if ~isiscalar(pad,'>-1'),
				error('Incorrect value for property ''pad'' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring_FMAT(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help MTSpectrum">MTSpectrum</a>'' for details).']);
	end
end

if isempty(frequency),
	if size(lfp,2) == 2,
		frequency = 1/median(diff(lfp(:,1)));
	else
		frequency = f;
	end
end

% Compute and plot spectrum
parameters.Fs = frequency;
if ~isempty(range), parameters.fpass = range; end
parameters.tapers = tapers;
parameters.pad = pad;
parameters.err = err;
[spectrum,f,s] = mtspectrumc(lfp(:,2),parameters);
s = s';
if strcmp(lower(show),'on'),
	PlotMean(f,log(spectrum),log(s(:,1)),log(s(:,2)),':');
end
