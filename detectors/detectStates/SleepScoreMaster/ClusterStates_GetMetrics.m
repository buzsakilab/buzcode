function [SleepScoreMetrics,StatePlotMaterials] = ClusterStates_GetMetrics(...
    basePath,SleepScoreLFP,EMG,overwrite,varargin)
%StateID(LFP,thLFP,EMG,sf_LFP,sf_EMG,figloc,WSEpisodes)
%   Detailed explanation goes here
%
%This function calucates the metrics for state scoring (EMG, Slow Wave
%Power, Theta). It uses SleepScoreLFP which is collected in
%PickSWTHChannel. And returns:
%   SleepScoreMetrics
%       .broadbandSlowWave
%       .thratio
%       .EMG
%       .t_clus     (time stamps, after ignored times have been removed)
%... also histograms and thresholds for subsequent scoring
%
%Dependencies: IDXtoINT
%
%Last Updated: 1/31/16
%DLevenstein
%% Params
p = inputParser;
addParameter(p,'onSticky',false)
addParameter(p,'ignoretime',[])
addParameter(p,'window',2)
addParameter(p,'smoothfact',15)
addParameter(p,'IRASA',true)
addParameter(p,'ThIRASA',true)
parse(p,varargin{:})
onSticky = p.Results.onSticky; 
ignoretime = p.Results.ignoretime; 
window = p.Results.window; 
smoothfact = p.Results.smoothfact; 
IRASA = p.Results.IRASA; 
ThIRASA = p.Results.ThIRASA; 

%This is the sticky trigger passed through to DetermineStates via histsandthreshs
if onSticky
    stickySW = true; stickyTH=false; stickyEMG=true;
else
    stickySW = false; stickyTH=false; stickyEMG=false;
end
%% Buzcode name of the SleepScoreMetrics.LFP.mat file
%Depreciated
 [datasetfolder,recordingname,extension] = fileparts(basePath);
% recordingname = [recordingname,extension]; % fileparts parses '.' into extension
% matfilename = fullfile(basePath,[recordingname,'.SleepScoreMetrics.LFP.mat']);
% plotmaterialsfilename = fullfile(basePath,[recordingname,'.StatePlotMaterials.mat']);
% 
% if exist(matfilename) & exist(plotmaterialsfilename) & overwrite == false
%     load(matfilename,'SleepScoreMetrics')
%     load(plotmaterialsfilename,'StatePlotMaterials')
%     return
% end


%Get the SW weights from SleepScoreLFP - if it's not there, use PSS
try
    SWweights = SleepScoreLFP.params.SWweights;
    SWfreqlist = SleepScoreLFP.params.SWfreqlist;
catch
    %Used to default to loading...  
    %load('SWweights.mat')
    SWweights = 'PSS';
end
%% Downsample and filter the LFP from PickSWTHChannel
%Make Downsample to niquest frequency

% if IRASA %IRASA needs sf>525
%     switch SleepScoreLFP.sf
%         case 1250 
%             downsamplefactor = 2;
%         otherwise
%             display('LFP is sampled too low, loading from .lfp file')
%             lfp = bz_GetLFP([SleepScoreLFP.SWchanID SleepScoreLFP.THchanID],...
%                  'basepath',basePath,'noPrompts',true,'downsample',2,...
%                  'intervals',SleepScoreLFP.t([1 end]));
%              swLFP = lfp.data(:,1);     thLFP = lfp.data(:,2);
%              t_LFP = lfp.timestamps;    sf_LFP = lfp.samplingRate;
%              downsamplefactor = 1;
%     end
% else
    switch SleepScoreLFP.sf
        case 1250
            downsamplefactor = 5;
        case 250
            downsamplefactor = 1;
        case 1000
            downsamplefactor = 4;
        otherwise
            display('sf not recognized... if only you made this able to set its own downsample...')
    end
% end

%Instead of taking SleepScoreLFP, this should use bz_getLFP and just needs
%SW/TH channels... can even get from sessionInfo.channeltags
swLFP = downsample(SleepScoreLFP.swLFP,downsamplefactor);
thLFP = downsample(SleepScoreLFP.thLFP,downsamplefactor);
t_LFP = downsample(SleepScoreLFP.t,downsamplefactor);
sf_LFP = SleepScoreLFP.sf/downsamplefactor;


%% Calculate broadbandslowwave metric

%Timing Parameters
noverlap = window-1; %1s dt

if strcmp(SWweights,'PSS')
    display('Calculating SW Mertric using Power Spectrum Slope')
    %Put the LFP in the right structure format
    lfp.data = swLFP;
    lfp.timestamps = t_LFP;
    lfp.samplingRate = sf_LFP;
    %Calculate PSS
    [specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,'frange',[4 90],'IRASA',IRASA);
    broadbandSlowWave = -specslope.data; %So NREM is higher as opposed to lower
    t_clus = specslope.timestamps;
    swFFTfreqs = specslope.freqs';
    specdt = 1./specslope.samplingRate;
    swFFTspec = 10.^spec.amp'; %To reverse log10 in bz_PowerSpectrumSlope

    IRASAsmooth = spec.IRASAsmooth';
    IRASAintercept = specslope.intercept;
    IRASAslope = specslope.data;
    
    % Remove transients before calculating SW histogram
    zFFTspec = NormToInt(log10(swFFTspec)','modZ');
    totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
    badtimes = find(totz>3);
else
    display(['Calculating SW Mertric using ',SleepScoreLFP.params.SWWeightsName])

    freqlist = logspace(0,2,100);
    [swFFTspec,swFFTfreqs,t_clus] = spectrogram(single(swLFP),window*sf_LFP,noverlap*sf_LFP,freqlist,sf_LFP);
    t_clus = t_clus'+t_LFP(1); %Offset for scoretime start
    swFFTspec = abs(swFFTspec);
    specdt = mode(diff(t_clus));
    zFFTspec = zscore(log10(swFFTspec)');
    % Remove transients before calculating SW histogram
    totz = zscore(abs(sum(zFFTspec')));
    badtimes = find(totz>5);

    %Calculate per-bin weights onto SlowWave
    assert(isequal(swFFTfreqs,SWfreqlist),...
        'spectrogram freqs.  are not what they should be...')
    broadbandSlowWave = zFFTspec*SWweights';
    
end

%Smooth and 0-1 normalize
broadbandSlowWave(badtimes) = nan;
broadbandSlowWave = smooth(broadbandSlowWave,smoothfact./specdt);
 
%% Calculate theta
%display('FFT Spectrum for Theta')

% %NarrowbandTheta
f_all = [2 20];
f_theta = [5 10];

if ThIRASA && strcmp(SWweights,'PSS')
    display('Calculating Theta Metric above PSS')
    %Put the LFP in the right structure format
    lfp.data = thLFP;
    lfp.timestamps = t_LFP;
    lfp.samplingRate = sf_LFP;
    %Calculate PSS
    [specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
        'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
    t_thclu = specslope.timestamps;
    specdt = 1./specslope.samplingRate;
    thFFTspec = specslope.resid';
    thFFTspec(thFFTspec<0)=0;
    
    IRASAsmooth_th = spec.IRASAsmooth';
    thFFTspec_raw = 10.^spec.amp';
    
    % Remove transients before calculating SW histogram
    zFFTspec = NormToInt(spec.amp,'modZ');
    totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
    badtimes_TH = find(totz>3);
    
    thFFTfreqs = specslope.freqs';
    thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
    thratio = max((thFFTspec(thfreqs,:)),[],1);
    
else
    display('Calculating Theta Metric using power ratio')
    freqlist = logspace(log10(f_all(1)),log10(f_all(2)),100);

    [thFFTspec,thFFTfreqs,t_thclu] = spectrogram(single(thLFP),window*sf_LFP,noverlap*sf_LFP,freqlist,sf_LFP);
    t_thclu = t_thclu+t_LFP(1); %Offset for scoretime start
    specdt = mode(diff(t_thclu));
    thFFTspec = (abs(thFFTspec));

    % Find transients for calculating TH
    zFFTspec = NormToInt(log10(thFFTspec)','modZ');
    totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
    badtimes_TH = find(totz>3);

    thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
    thpower = sum((thFFTspec(thfreqs,:)),1);
    allpower = sum((thFFTspec),1);

    thratio = thpower./allpower;    %Narrowband Theta
end

thratio(badtimes_TH) = nan;     %Remove transients
thratio = smooth(thratio,smoothfact./specdt);


%% Remove ignoretimes (after smoothing), before normalizoing
if ~isempty(ignoretime) || ~isempty(badtimes_TH) || ~isempty(badtimes)
	ignoretimeIDX = InIntervals(t_clus,ignoretime) | isnan(broadbandSlowWave) | isnan(thratio);
    broadbandSlowWave(ignoretimeIDX) = [];
    thratio(ignoretimeIDX) = [];
    t_clus(ignoretimeIDX) = [];
    t_thclu(ignoretimeIDX) = [];
end
broadbandSlowWave = bz_NormToRange(broadbandSlowWave,[0 1]);
thratio = bz_NormToRange(thratio,[0 1]);
 
%% EMG
dtEMG = 1/EMG.samplingFrequency;
EMG.smoothed = smooth(EMG.data,smoothfact/dtEMG,'moving');

%remove any t_clus before/after t_emg
prEMGtime = t_clus<EMG.timestamps(1) | t_clus>EMG.timestamps(end);
broadbandSlowWave(prEMGtime) = []; thratio(prEMGtime) = []; t_clus(prEMGtime) = [];

%interpolate to FFT time points;
EMG = interp1(EMG.timestamps,EMG.smoothed,t_clus,'nearest');

%Min/Max Normalize
EMG = bz_NormToRange(EMG,[0 1]);


%% BELOW IS GETTING HISTOGRAMS AND THRESHOLDS FOR SCORING IN         %%
%  CLUSTERSTATES_DetermineStates and visualization in TheStateEditor %%\

%% Divide PC1 for SWS
%DL Note: should replace all of this with calls to bz_BimodalThresh
numpeaks = 1;
numbins = 12;
%numbins = 12; %for Poster...
while numpeaks ~=2
    [swhist,swhistbins]= hist(broadbandSlowWave,numbins);
    
    [PKS,LOCS] = findpeaks_SleepScore(swhist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

SWdiptest = bz_hartigansdiptest(broadbandSlowWave);
betweenpeaks = swhistbins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks_SleepScore(-swhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

swthresh = betweenpeaks(diploc);


 
%SWS time points
NREMtimes = (broadbandSlowWave >swthresh);


%% Then Divide EMG
numpeaks = 1;
numbins = 12;
if sum(isnan(EMG))>0
   error('EMG seems to have NaN values...') 
end

while numpeaks ~=2
    %[EMGhist,EMGhistbins]= hist(EMG(NREMtimes==0),numbins); 
    [EMGhist,EMGhistbins]= hist(EMG,numbins); %changed back 6/18/19

    [PKS,LOCS] = findpeaks_SleepScore([0 EMGhist],'NPeaks',2);
    LOCS = sort(LOCS)-1;
    numbins = numbins+1;
    numpeaks = length(LOCS);
    
    if numpeaks ==100
        display('Something is wrong with your EMG')
        return
    end
end

EMGdiptest = bz_hartigansdiptest(EMG);
betweenpeaks = EMGhistbins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks_SleepScore(-EMGhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

EMGthresh = betweenpeaks(diploc);

MOVtimes = (broadbandSlowWave(:)<swthresh & EMG(:)>EMGthresh);


%% Then Divide Theta (During NonMoving)
numpeaks = 1;
numbins = 12;
while numpeaks ~=2 && numbins <=25
    %[THhist,THhistbins]= hist(thratio(SWStimes==0 & MOVtimes==0),numbins);
    [THhist,THhistbins]= hist(thratio(MOVtimes==0),numbins);

    [PKS,LOCS] = findpeaks_SleepScore(THhist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

THdiptest = bz_hartigansdiptest(thratio(MOVtimes==0));

if numpeaks ~= 2
	display('No bimodal dip found in theta. Trying to exclude NREM...')

    numbins = 12;
    %numbins = 15; %for Poster...
    while numpeaks ~=2 && numbins <=25
        [THhist,THhistbins]= hist(thratio(NREMtimes==0 & MOVtimes==0),numbins);

        [PKS,LOCS] = findpeaks_SleepScore(THhist,'NPeaks',2,'SortStr','descend');
        LOCS = sort(LOCS);
        numbins = numbins+1;
        numpeaks = length(LOCS);

    end
    
	try
        THdiptest = bz_hartigansdiptest(thratio(NREMtimes==0 & MOVtimes==0));
    catch
	end
end

if length(PKS)==2
    betweenpeaks = THhistbins(LOCS(1):LOCS(2));
    [dip,diploc] = findpeaks_SleepScore(-THhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

    THthresh = betweenpeaks(diploc);

    REMtimes = (broadbandSlowWave<swthresh & EMG<EMGthresh & thratio>THthresh);
else
    display('No bimodal dip found in theta. Use TheStateEditor to manually select your threshold (hotkey: A)')
    THthresh = 0;
%     REMtimes =(broadbandSlowWave<swthresh & EMG<EMGthresh);
end

histsandthreshs = v2struct(swhist,swhistbins,swthresh,EMGhist,EMGhistbins,...
    EMGthresh,THhist,THhistbins,THthresh,...
    stickySW,stickyTH,stickyEMG);

%% Ouput Structure: StateScoreMetrics
LFPparams = SleepScoreLFP.params;
WindowParams.window = window;
WindowParams.smoothwin = smoothfact;
THchanID = SleepScoreLFP.THchanID; SWchanID = SleepScoreLFP.SWchanID;

SleepScoreMetrics = v2struct(broadbandSlowWave,thratio,EMG,...
    t_clus,badtimes,badtimes_TH,histsandthreshs,LFPparams,WindowParams,THchanID,SWchanID,...
    recordingname,THdiptest,EMGdiptest,SWdiptest);
%save(matfilename,'SleepScoreMetrics');

StatePlotMaterials = v2struct(t_clus,swFFTfreqs,swFFTspec,thFFTfreqs,thFFTspec);
if exist('IRASAsmooth_th','var')
    StatePlotMaterials.IRASAsmooth_th = IRASAsmooth_th;
    StatePlotMaterials.thFFTspec_raw = thFFTspec_raw;

end
if exist('IRASAsmooth','var')
    StatePlotMaterials.IRASAsmooth = IRASAsmooth;
    StatePlotMaterials.IRASAintercept = IRASAintercept;
    StatePlotMaterials.IRASAslope = IRASAslope;
end
%save(plotmaterialsfilename,'StatePlotMaterials'); 
