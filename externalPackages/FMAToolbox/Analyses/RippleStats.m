function [maps,data,stats] = RippleStats(filtered,ripples,varargin)

%RippleStats - Compute descriptive stats for ripples (100~200Hz oscillations).
%
%  USAGE
%
%    [maps,data,stats] = RippleStats(filtered,ripples,<options>)
%
%    filtered       ripple-band filtered samples (one channel)
%    ripples        ripple timing information (obtained using <a href="matlab:help FindRipples">FindRipples</a>)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   sampling rate (in Hz) (default = 1250Hz)
%     'durations'   durations before and after ripple peak (in s)
%                   (default = [-0.5 0.5])
%    =========================================================================
%
%  OUTPUT
%
%    maps.ripples               instantaneous amplitude (one ripple per row)
%    maps.frequency             instantaneous frequency
%    maps.phase                 instantaneous phase
%    maps.amplitude             enveloppe amplitude
%    data.peakFrequency         frequency at peak
%    data.peakAmplitude         amplitude at peak
%    data.duration              durations
%    stats.durationAmplitude    correlation between duration and amplitude (rho, p)
%    stats.durationFrequency    correlation between duration and frequency (rho, p)
%    stats.amplitudeFrequency   correlation between amplitude and frequency (rho, p)
%    stats.acg.data             autocorrelogram data
%    stats.acg.t                autocorrelogram time bins
%
%  SEE
%
%    See also FindRipples, SaveRippleEvents, PlotRippleStats.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
samplingRate = 1250;
durations = [-50 50]/1000;
nCorrBins = 50;
%  corrBinSize = 400;
corrDuration = 20;
corrBinSize = 0.01;

% Check number of parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help RippleStats">RippleStats</a>'' for details).');
end

% Check parameter sizes
if size(filtered,2) ~= 2,
	error('Parameter ''filtered'' is not a Nx2 matrix (type ''help <a href="matlab:help RippleStats">RippleStats</a>'' for details).');
end
%  if size(ripples,2) ~= 4,
%  	error('Parameter ''ripples'' is not a Nx4 matrix (type ''help <a href="matlab:help RippleStats">RippleStats</a>'' for details).');
%  end
if size(ripples,2) < 3,
	error('Parameter ''ripples'' is not a Nx4 or Nx4 matrix (type ''help <a href="matlab:help RippleStats">RippleStats</a>'' for details).');
end



% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help RippleStats">RippleStats</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'frequency',
			samplingRate = varargin{i+1};
			if ~isdscalar(samplingRate,'>0'),
				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help RippleStats">RippleStats</a>'' for details).');
			end
		case 'durations',
			durations = varargin{i+1};
			if ~isdvector(durations,'#2','<'),
				error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help RippleStats">RippleStats</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help RippleStats">RippleStats</a>'' for details).']);
	end
end

nBins = floor(samplingRate*diff(durations)/2)*2+1; % must be odd
nHalfCenterBins = 3;
centerBin = ceil(nBins/2);
centerBins = centerBin-nHalfCenterBins:centerBin+nHalfCenterBins;

% Compute instantaneous phase and amplitude
h = hilbert(filtered(:,2));
phase(:,1) = filtered(:,1);
phase(:,2) = angle(h);
amplitude(:,1) = filtered(:,1);
amplitude(:,2) = abs(h);
unwrapped(:,1) = filtered(:,1);
unwrapped(:,2) = unwrap(phase(:,2));
% Compute instantaneous frequency
frequency = Diff(unwrapped,'smooth',0);
frequency(:,2) = frequency(:,2)/(2*pi);

% Compute ripple map
[r,i] = Sync(filtered,ripples(:,2),'durations',durations);
maps.ripples = SyncMap(r,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute frequency Map
[f,i] = Sync(frequency,ripples(:,2),'durations',durations);
maps.frequency = SyncMap(f,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute phase map
[p,i] = Sync(phase,ripples(:,2),'durations',durations);
maps.phase = SyncMap(p,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute amplitude map
[a,i] = Sync(amplitude,ripples(:,2),'durations',durations);
maps.amplitude = SyncMap(a,i,'durations',durations,'nbins',nBins,'smooth',0);

% Ripple frequency and amplitude at peak
data.peakFrequency = maps.frequency(:,centerBin);
data.peakAmplitude = maps.amplitude(:,centerBin);

% Ripple durations
data.duration = ripples(:,3)-ripples(:,1);

% Autocorrelogram and correlations
%  if nargin > 2,
%  	[stats.acg.data,stats.acg.t] = CCG(ripples(:,2),1,'binSize',corrBinSize,'halfBins',nCorrBins/2);
	[stats.acg.data,stats.acg.t] = CCG(ripples(:,2),1,'binSize',corrBinSize);
	[stats.amplitudeFrequency.rho,stats.amplitudeFrequency.p] = corrcoef(data.peakAmplitude,data.peakFrequency);
	[stats.durationFrequency.rho,stats.durationFrequency.p] = corrcoef(data.duration,data.peakFrequency);
	[stats.durationAmplitude.rho,stats.durationAmplitude.p] = corrcoef(data.duration,data.peakAmplitude);
%  end
