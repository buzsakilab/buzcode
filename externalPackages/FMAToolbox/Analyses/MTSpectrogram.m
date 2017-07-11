function [spectrogram,t,f] = MTSpectrogram(lfp,varargin)

%MTSpectrogram - Compute LFP spectrogram by multi-taper estimation.
%
%  The spectrogram is computed using the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  USAGE
%
%    [spectrogram,t,f] = MTSpectrogram(lfp,<options>)
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
%     'window'      duration (in s) of the time window (default = 5)
%     'overlap'     overlap between successive windows (default = window/2)
%     'step'        step between successive windows (default = window/2)
%     'tapers'      relative resolution and order of the tapers [NW K]
%                   (default = [3 5])
%     'pad'         FFT padding (see help for <a href="matlab:help mtspecgramc">mtspecgramc</a>) (default = 0)
%     'show'        plot spectrogram (default = 'off')
%     'cutoffs'     cutoff values for color plot (default = [0 13])
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
%    spectrogram    time-frequency matrix
%    t              time bins
%    f              frequency bins
%
%  DEPENDENCIES
%
%    This function requires the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  SEE
%
%    See also MTSpectrum, SpectrogramBands, MTCoherence, MTCoherogram, PlotColorMap.

% Copyright (C) 2004-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure chronux is installed and functional
CheckChronux();

% Defaults
f = 1250;
frequency = [];
window = 5;
range = [];
overlap = [];
step = [];
show = 'off';
cutoffs = [0 13];
tapers = [3 5];
pad = 0;

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
end

% Check parameter sizes
if size(lfp,2) ~= 1 && size(lfp,2) ~= 2,
	error('Parameter ''lfp'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'frequency',
			frequency = varargin{i+1};
			if ~isdscalar(frequency,'>0'),
				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'range',
			range = varargin{i+1};
			if ~isdvector(range,'#2','>=0','<'),
				error('Incorrect value for property ''range'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'window',
			window = varargin{i+1};
			if ~isdscalar(window,'>0'),
				error('Incorrect value for property ''window'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'overlap',
			overlap = varargin{i+1};
			if ~isdscalar(overlap,'>0'),
				error('Incorrect value for property ''overlap'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'step',
			step = varargin{i+1};
			if ~isdscalar(step,'>0'),
				error('Incorrect value for property ''step'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'tapers',
			tapers = varargin{i+1};
			if ~isivector(tapers,'#2','>0'),
				error('Incorrect value for property ''tapers'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'pad',
			pad = varargin{i+1};
			if ~isiscalar(pad,'>-1'),
				error('Incorrect value for property ''pad'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring_FMAT(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'cutoffs',
			cutoffs = varargin{i+1};
			if ~isdvector(cutoffs,'#2','>=0','<'),
				error('Incorrect value for property ''cutoffs'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).']);
	end
end

% Determine LFP frequency
if isempty(frequency),
	if size(lfp,2) == 2,
		frequency = 1/median(diff(lfp(:,1)));
	else
		frequency = f;
	end
end

% Determine step/overlap
if isempty(step),
	if isempty(overlap),
		overlap = window/2;
	end
else
	if isempty(overlap),
		overlap = window-step;
	elseif overlap ~= window-step,
		error('Incompatible ''step'' and ''overlap'' parameters (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
	end
end

% If lfp was given without timestamps, add them...
if size(lfp,2) == 1
    lfp(:,2) = lfp(:,1);
    lfp(:,1) = [1:size(lfp,1)]./frequency;
end

% Compute and plot spectrogram
parameters.Fs = frequency;
if ~isempty(range), parameters.fpass = range; end
parameters.tapers = tapers;
parameters.pad = pad;
[spectrogram,t,f] = mtspecgramc(lfp(:,2),[window window-overlap],parameters);
t = t'+lfp(1,1);
f = f';
spectrogram = spectrogram';
if strcmp(lower(show),'on'),
%  	figure;
	logTransformed = log(spectrogram);
	PlotColorMap(logTransformed,1,'x',t,'y',f,'cutoffs',cutoffs,'newfig','off');
	xlabel('Time (s)');
	ylabel('Frequency (Hz)');
	title('Power Spectrogram');
end
