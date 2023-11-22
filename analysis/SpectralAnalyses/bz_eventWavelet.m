
function [ wavAvg, lfpAvg ] = bz_eventWavelet (lfp, events, varargin)

% [ wavAvg, lfpAvg ] = bz_eventWavelet (lfp, events, varargin)
% Calculates event-triggered (i.e. SWRs) wavelet spectrogram

% INPUT
%    lfp            a buzcode structure with fields lfp.data (only 1 channel)
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   -lfp can also be a [t x 1] timeseries signal. in which
%                   case you need to input 'samplingRate'
%   events          events timestamps (in sec)

%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       twin        time window around events to calculate average. Default [0.1 0.1]
%       plotWave    true/false. Default true.
%       plotLFP     true/false. Default true.
%       whitening   perform lfp whitening. Default true.
%       tsmth       degree of temporal smoothing in s. Default 0.0001 (0=no smothing).
%       waveDb      spectrum in decibel units. Default false.

%   To pass them to bz_WaveSpec:
%       'frange'    [low frequency, high frequency]     (default: [50 200])
%       'nfreqs'    number of frequencies               (default: 100
%       'roundfreqs' round freqs to unique integer vals (default: false
%                       *Note this may decrease number of freqs
%       'nfreqs'    number of frequencies               (default: 100
%       'ncyc'      number of cycles in your wavelet    (default: 5)
%       'fvector'   predefined vector of frequencies 
%       'space'     'log' or 'lin'  spacing of f's      (default: 'log')
%       'samplingRate' (only if input is not a buzcode structure)

%    =========================================================================

% OUTPUT:
%    wavAvg           a buzcode structure with fields wavAvg.data,
%                                                   wavAvg.timestamps
%                                                   wavAvg.samplingRate
%                                                   wavAvg.params
%    lfpAvg [        a buzcode structure with fields lpf.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                                                   lfp.channels 
%                                                   lfp.params

% Antonio FR, 5/20

% TODO:
%    - z-scoring of power based on session baseline
%    - allow reading from binary lfp file instead of from a lfp.mat
%    - improve plot 

%% Parse inputs

p = inputParser;
addParameter(p,'channels',1:size(lfp.data,2),@isvector);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'twin',[0.1 0.1],@isnumeric);
addParameter(p,'plotWave',true,@islogical);
addParameter(p,'plotLFP',true,@islogical);
addParameter(p,'whitening',true,@islogical);
addParameter(p,'tsmth',0.0001,@isnumeric);
addParameter(p,'waveDb',false,@islogical);
addParameter(p,'frange',[50 250],@isnumeric);
addParameter(p,'nfreqs',100,@isnumeric);
addParameter(p,'ncyc',5,@isnumeric);
addParameter(p,'space','log');
addParameter(p,'roundfreqs',false,@islogical);
addParameter(p,'fvector',[]);

parse(p,varargin{:});
channels = p.Results.channels;
samplingRate = p.Results.samplingRate;
plotWave = p.Results.plotWave;
plotLFP = p.Results.plotLFP;
whitening = p.Results.whitening;
tsmth = p.Results.tsmth;
waveDb = p.Results.waveDb;
frange = p.Results.frange;
nfreqs = p.Results.nfreqs;
ncyc = p.Results.ncyc;
space = p.Results.space;
roundfreqs = p.Results.roundfreqs;
fvector = p.Results.fvector;

%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
end

twin = p.Results.twin*samplingRate;
events = round(events*samplingRate);

%% Conpute event-triggered LFP average

events = events((events + twin(2) <= size(data,1)) & (events - twin(1) > 0));
lfp_temp = nan(twin(1)+twin(2)+1,1,length(events));

for e = 1:length(events)
    lfp_temp(:,:,e) = data(events(e)-twin(1):events(e)+twin(2),1);
end

lfp_avg = nanmean(lfp_temp,3);

%% Whiten LFP

if whitening
%     temp = WhitenSignal( lfp.data', [], 0 )'; %  NEED TO MAKE BUZCODE FUNCTION
%     lfp.data = int16(temp'); clear temp;
    lfpwhiten = bz_whitenLFP(lfp);
    lfp.data = lfpwhiten.data; clear temp;    
end

%% temporal smoothing

if tsmth > 0
    winF = gausskernel(0, ceil(tsmth *samplingRate),1, 6* ceil(tsmth *samplingRate) + 1 );
    winF = winF / sum( winF(:) );

    if length( winF ) > 1
        C = length(winF);
        D = ceil(C/2) - 1;
        temp = filter(winF,1,[flipud(lfp.data(1:C,:));lfp.data;flipud(lfp.data(end-C+1:end,:))]);
        temp = temp(1+C+D:end-C+D,:);
        lfp.data = temp; clear temp;
    end
end

%% Conpute wavelets

for i = 1:length(events)
    eventWin(i,1) = lfp.timestamps(events(i)-twin(1));
    eventWin(i,2) = lfp.timestamps(events(i)+twin(2));
end

wavespec = bz_WaveSpec(lfp,'intervals',eventWin,'frange',frange,'nfreqs',nfreqs,'roundfreqs',roundfreqs,...
            'nfreqs',nfreqs,'ncyc',ncyc,'fvector',fvector,'space',space);

% Pull out individual wavelet spectrograms
win = twin(1)+twin(2)+1;
wavespec_ints = nan( win, nfreqs,  length(eventWin) );
for kp = 1:size(eventWin,1)
    tt = InIntervals( wavespec.timestamps, eventWin(kp,:) );
    wavespec_ints(1:sum(tt),:,kp) = wavespec.data( tt , :);
end

% Power spectrum 
wavespec_ints_amp = abs(wavespec_ints);
if waveDb
   avgSpect =  20*log10(squeeze(nanmean(wavespec_ints_amp,3)));
else
   avgSpect =  squeeze(nanmean(wavespec_ints_amp,3));
end

% Peak freq
[t,f] = find(max(avgSpect(:))==(avgSpect));
wavAvg.maxFreq = wavespec.freqs(f);
wavAvg.maxPow = max(avgSpect(:));

%% generate output structure

wavAvg.data = avgSpect;
wavAvg.timestamps = -p.Results.twin(1):(1/samplingRate):p.Results.twin(2);
wavAvg.freqs = wavespec.freqs;
wavAvg.nfreqs = wavespec.nfreqs;
wavAvg.samplingRate = wavespec.samplingRate;
if isfield(wavespec,'channels')
    wavAvg.channels = wavespec.channels;
end
wavAvg.params = wavespec.filterparms;
wavAvg.params.whitening = whitening;
wavAvg.params.tsmth = tsmth;

lfpAvg.data = lfp_avg;
lfpAvg.timestamps = -twin(1):twin(2);
lfpAvg.samplingRate = samplingRate;
lfpAvg.channels = channels; 

%% Plot
if plotWave
   figure;
   contourf(wavAvg.timestamps*1000,wavAvg.freqs,wavAvg.data',30,'LineColor','none');hold on;
   set(gca,'YScale','log');
   ylim([wavespec.freqs(1) wavespec.freqs(end)]);
   colormap jet;
   xlabel('time (ms)'); ylabel('frequency (Hz)');
    
end
  
if plotLFP
   fmin =  wavespec.freqs(round(length(wavespec.freqs)/4)); 
   fmax =  fmin*2;
   lfpSc =(lfp_avg-min(lfp_avg))*(fmax-fmin)/(max(lfp_avg)-min(lfp_avg))+fmin;
   plot(wavAvg.timestamps*1000,lfpSc,'w','LineWidth',2);hold on;
    
end

end
