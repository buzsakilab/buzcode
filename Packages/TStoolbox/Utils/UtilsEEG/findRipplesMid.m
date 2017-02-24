function midRip = findRipplesMid(eeg,epoch)

% findRipplesMid finds ripples peaks in hippocampal signals.
% USAGE:
%  midRip = findRipplesMid(eeg,epoch)
%  input:
%  	eeg: tsd of unfiltered eeg signal
%  	epoch: intervalset on which ripples have to be searched (typically: SWS epochs). 
%  output:
%  	midRip: a tsd of ripples peaks times and ripples peaks amplitude (in z-score of the filtered signal)
%  
%  Important: don't restrict eeg signals to SWS before as discontinuities in signal will yield detection of erroneous ripples (e.g. noise)
%
%  Adrien Peyrache, 2007

%Parameters:
thresRip = 2; % threshold for enveloppe detection
minPeaks = 6; % threshold for peaks detection (in z-score)

dt = median(diff(Range(eeg,'s')));
Fn = 1/(2*dt);
b = fir1(96,[100 300]/Fn);

dEegFilt = filtfilt(b,1,Data(eeg));
eegFilt = tsd(Range(eeg),dEegFilt);
eegFilt = Restrict(eegFilt,epoch);;
dEegFilt = Data(eegFilt);
stdEeg = std(dEegFilt);
dEegFilt = dEegFilt/stdEeg;
eegFilt = tsd(Range(eegFilt),dEegFilt);

dEegFAbs = abs(dEegFilt);

gw = gausswin(8*Fn/200);
gw = gw/sum(gw);

dEegFAbSm = convn(dEegFAbs,gw,'same');
eegSm = tsd(Range(eegFilt),dEegFAbSm);

ripInt = thresholdIntervals(eegSm,thresRip,'Direction','Above');
ripInt = mergeCloseIntervals(ripInt,200); %Old 400
ripInt = dropShortIntervals(ripInt,200);
ripInt = dropLongIntervals(ripInt,2000); %old 1000
nbRipples = length(Start(ripInt)),

maxRip = max(eegFilt,ripInt);
maxRipTime = Range(maxRip);
maxRipAmp = Data(maxRip);
ix = maxRipAmp>minPeaks;
maxRipTime = maxRipTime(ix);
maxRipAmp = maxRipAmp(ix);

[dummy,sortIx] = sort(maxRipTime);
midRip = tsd(maxRipTime(sortIx),maxRipAmp(sortIx));
