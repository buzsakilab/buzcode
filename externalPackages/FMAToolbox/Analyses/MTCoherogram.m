function [coherogram,phase,t,f] = MTCoherogram(lfp1,lfp2,varargin)

%MTCoherogram - Compute LFP coherogram by multi-taper estimation.
%
%  USAGE
%
%    [coherogram,phase,t,f] = MTCoherogram(lfp1,lfp2,<options>)
%
%    lfp1,lfp2      wide-band LFPs (one channel each).
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
%     'pad'         FFT padding (see help for <a href="matlab:help cohgramc">cohgramc</a>) (default = 0)
%     'show'        plot results (default = 'off')
%     'cutoffs'     cutoff values for color plot (default = [0 1])
%    =========================================================================
%
%  NOTES
%
%    The LFP can be provided either as a time stamped matrix (list of time-voltage
%    pairs), or as a voltage vector - in which case the frequency must be specified.
%
%    The time displacement between successive short time coherences can be supplied
%    either as a 'step' (explicit time difference) or as an 'overlap' (between
%    successive time windows).
%
%  OUTPUT
%
%    coherogram     coherogram magnitude
%    phase          coherogram phase
%    t              time bins
%    f              frequency bins
%
%  DEPENDENCIES
%
%    This function requires the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  SEE
%
%    See also MTCoherence, MTSpectrum, MTSpectrogram, PlotColorMap.

% Copyright (C) 2010-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
warning('This function has been deprecated, and will be removed in the future. Try using bz_MTCoherogram.m')
% Make sure chronux is installed and functional
CheckChronux('cohgramc');

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
cutoffs = [0 1];

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
end

% Check parameter sizes
if size(lfp1,2) ~= 1 && size(lfp1,2) ~= 2,
	error('Parameter ''lfp1'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
end
if size(lfp2,2) ~= 1 && size(lfp2,2) ~= 2,
	error('Parameter ''lfp2'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'frequency',
			frequency = varargin{i+1};
			if ~isdscalar(frequency,'>0'),
				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
			end
		case 'range',
			range = varargin{i+1};
			if ~isdvector(range,'#2','<','>=0'),
				error('Incorrect value for property ''range'' (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
			end
		case 'window',
			window = varargin{i+1};
			if ~isdscalar(window,'>0'),
				error('Incorrect value for property ''window'' (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
			end
		case 'overlap',
			overlap = varargin{i+1};
			if ~isdscalar(overlap,'>0'),
				error('Incorrect value for property ''overlap'' (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
			end
		case 'step',
			step = varargin{i+1};
			if ~isdscalar(step,'>0'),
				error('Incorrect value for property ''step'' (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
			end
		case 'tapers',
			tapers = varargin{i+1};
			if ~isivector(tapers,'#2','>0'),
				error('Incorrect value for property ''tapers'' (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
			end
		case 'pad',
			pad = varargin{i+1};
			if ~isiscalar(pad,'>-1'),
				error('Incorrect value for property ''pad'' (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring_FMAT(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
			end
		case 'cutoffs',
			cutoffs = varargin{i+1};
			if ~isdvector(cutoffs,'#2','>=0','<'),
				error('Incorrect value for property ''cutoffs'' (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).']);
	end
end

% Determine LFP frequency
if isempty(frequency),
	if size(lfp1,2) == 2,
		frequency = 1/median(diff(lfp1(:,1)));
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
		error('Incompatible ''step'' and ''overlap'' parameters (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
	end
end

% Compute and plot coherogram
parameters.Fs = frequency;
if ~isempty(range), parameters.fpass = range; end
parameters.tapers = tapers;
parameters.pad = pad;
[coherogram,phase,~,~,~,t,f] = cohgramc(lfp1(:,2),lfp2(:,2),[window window-overlap],parameters);
t = t'+lfp1(1,1);
f = f';
coherogram = coherogram';
% coherogram = permute(coherogram,[2 1 3]);  % Previous code by Gabrielle Girardeau, keep it around just in case
phase = phase';
if strcmp(lower(show),'on'),
  figure;hold on;
  subplot(2,1,1);
  PlotColorMap(coherogram,'x',t,'y',f,'cutoffs',cutoffs,'newfig','off');
  xlabel('Time (s)');
  ylabel('Frequency (Hz)');
  title('Coherogram Amplitude');
  subplot(2,1,2);
  PlotColorMap(phase,'x',t,'y',f,'cutoffs',[-pi pi],'newfig','off');
  xlabel('Time (s)');
  ylabel('Frequency (Hz)');
  title('Coherogram Phase');
end
