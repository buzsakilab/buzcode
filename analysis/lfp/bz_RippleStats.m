function [maps,data,stats] = bz_RippleStats(filtered,timestamps,ripples,varargin)

%RippleStats - Compute descriptive stats for ripples (100~200Hz oscillations).
%
%  USAGE
%
%    [maps,data,stats] = bz_RippleStats(filtered,timestamps,ripples,<options>)
%
%    filtered       ripple-band filtered samples (one channel)
%    timestamps	    timestamps (in seconds) to match filtered vector
%    ripples        ripple timing information (STRUCT VERSION) (obtained using <a href="matlab:help bz_FindRipples">bz_FindRipples</a>)
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

% edited by David Tingley to fit buzcode formatting standards, 2017

% Default values
samplingRate = 1250;
durations = [-75 75]/1000;
nCorrBins = 50;
%  corrBinSize = 400;
corrDuration = 20;
corrBinSize = 0.01;

% Check number of parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help RippleStats">RippleStats</a>'' for details).');
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

idx = ceil((ripples.peaks-ripples.timestamps(:,1))*ripples.detectorinfo.detectionparms.frequency);

% Compute instantaneous phase and amplitude
h = hilbert(filtered);
phase = angle(h);
amplitude = abs(h);
unwrapped = unwrap(phase);
% Compute instantaneous frequency
frequency = bz_Diff(medfilt1(unwrapped,12),timestamps,'smooth',0);
frequency = frequency/(2*pi);

% Compute ripple map
[r,i] = Sync([timestamps filtered],ripples.peaks,'durations',durations);
maps.ripples = SyncMap(r,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute frequency Map
[f,i] = Sync([timestamps frequency],ripples.peaks,'durations',durations);
maps.frequency = SyncMap(f,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute phase map
[p,i] = Sync([timestamps phase],ripples.peaks,'durations',durations);
maps.phase = SyncMap(p,i,'durations',durations,'nbins',nBins,'smooth',0);

% Compute amplitude map
[a,i] = Sync([timestamps amplitude],ripples.peaks,'durations',durations);
maps.amplitude = SyncMap(a,i,'durations',durations,'nbins',nBins,'smooth',0);

idx(idx>length(maps.frequency(1,:))) = length(maps.frequency(1,:));
% Ripple frequency and amplitude at peak
% for i= 1:length(ripples.timestamps)
%     data.peakFrequency(i) = maps.frequency(i,idx(i));
%     data.peakAmplitude(i) = maps.amplitude(i,idx(i));
% end
% data.peakFrequency = data.peakFrequency';
% data.peakAmplitude = data.peakAmplitude';

data.peakFrequency = maps.frequency(:,centerBin);
data.peakAmplitude = maps.amplitude(:,centerBin);

% Ripple durations
data.duration = abs(diff(ripples.timestamps'))';

% Autocorrelogram and correlations
%  if nargin > 2,
%  	[stats.acg.data,stats.acg.t] = CCG(ripples.peaks,1,'binSize',corrBinSize,'halfBins',nCorrBins/2);
	[stats.acg.data,stats.acg.t] = CCG(ripples.peaks,ones(length(ripples.peaks),1),'binSize',corrBinSize);
	[stats.amplitudeFrequency.rho,stats.amplitudeFrequency.p] = corrcoef(data.peakAmplitude,data.peakFrequency);
	[stats.durationFrequency.rho,stats.durationFrequency.p] = corrcoef(data.duration,data.peakFrequency);
	[stats.durationAmplitude.rho,stats.durationAmplitude.p] = corrcoef(data.duration,data.peakAmplitude);
%  end
