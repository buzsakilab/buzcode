function [SleepScoreMetrics,StatePlotMaterials] = ClusterStates_GetMetrics(...
    basePath,SleepScoreLFP,EMG,overwrite)
%StateID(LFP,thLFP,EMG,sf_LFP,sf_EMG,figloc,WSEpisodes)
%   Detailed explanation goes here
%
%
%
%Dependencies: IDXtoINT
%
%Last Updated: 1/31/16
%DLevenstein

%% Buzcode name of the SleepScoreMetrics.LFP.mat file
[datasetfolder,recordingname,extension] = fileparts(basePath);
recordingname = [recordingname,extension]; % fileparts parses '.' into extension
matfilename = fullfile(basePath,[recordingname,'.SleepScoreMetrics.LFP.mat']);
plotmaterialsfilename = fullfile(basePath,[recordingname,'.StatePlotMaterials.mat']);

if exist(matfilename) & exist(plotmaterialsfilename) & overwrite == false
    load(matfilename,'SleepScoreMetrics')
    load(plotmaterialsfilename,'StatePlotMaterials')
    return
end

%% Downsample and filter
%Make Downsample to niquest frequency

if SleepScoreLFP.sf == 1250
    downsamplefactor = 5;
elseif SleepScoreLFP.sf == 250
    downsamplefactor = 1;
elseif SleepScoreLFP.sf == 1000
    downsamplefactor = 4;
else
    display('sf not recognized... if only you made this able to set its own downsample...')
end
swLFP = downsample(SleepScoreLFP.swLFP,downsamplefactor);
thLFP = downsample(SleepScoreLFP.thLFP,downsamplefactor);
sf_LFP = SleepScoreLFP.sf/downsamplefactor;


%% Calculate Spectrogram
%display('FFT Spectrum for Broadband LFP')

freqlist = logspace(0,2,100);
window = 10;   %s
noverlap = 9;  %s
window = window*sf_LFP;
noverlap = noverlap*sf_LFP;
[swFFTspec,swFFTfreqs,t_clus] = spectrogram(single(swLFP),window,noverlap,freqlist,sf_LFP);
swFFTspec = abs(swFFTspec);
[zFFTspec,mu,sig] = zscore(log10(swFFTspec)');

%% Remove transients before calculating SW histogram
%this should be it's own whole section - removing/detecting transients
totz = zscore(abs(sum(zFFTspec')));
badtimes = find(totz>5);
zFFTspec(badtimes,:) = 0;
 
%% Set Broadband filter weights for Slow Wave
load('SWweights.mat')
assert(isequal(freqlist,SWfreqlist), 'spectrogram freqs.  are not what they should be...')
broadbandSlowWave = zFFTspec*SWweights';
 
%% Smooth and 0-1 normalize
smoothfact = 10; %units of si_FFT
thsmoothfact = 10; %used to be 15

broadbandSlowWave = smooth(broadbandSlowWave,smoothfact);
broadbandSlowWave = (broadbandSlowWave-min(broadbandSlowWave))./max(broadbandSlowWave-min(broadbandSlowWave));

 
%% Calculate theta
%display('FFT Spectrum for Theta')

% %NarrowbandTheta
f_all = [2 20];
f_theta = [5 10];
freqlist = logspace(log10(f_all(1)),log10(f_all(2)),100);


[thFFTspec,thFFTfreqs] = spectrogram(single(thLFP),window,noverlap,freqlist,sf_LFP);
thFFTspec = (abs(thFFTspec));
[~,mu_th,sig_th] = zscore(log10(thFFTspec)');

thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
allpower = sum((thFFTspec),1);
thpower = sum((thFFTspec(thfreqs,:)),1);

thratio = thpower./allpower;    %Narrowband Theta
thratio = smooth(thratio,thsmoothfact);
thratio = (thratio-min(thratio))./max(thratio-min(thratio));
 
%% EMG
dtEMG = 1/EMG.samplingFrequency;
t_EMG = (1:length(EMG.data))*dtEMG;
EMG = smooth(EMG.data,smoothfact/dtEMG);
EMG = (EMG-min(EMG))./max(EMG-min(EMG));

reclength = round(t_EMG(end));

%downsample to FFT time points;
[~,t_intersect] = intersect(t_EMG,t_clus);
EMG = EMG(t_intersect);
t_EMG = t_EMG(t_intersect);



%% Divide PC1 for SWS
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


betweenpeaks = swhistbins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks_SleepScore(-swhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

swthresh = betweenpeaks(diploc);

%Set transients to wake state
broadbandSlowWave(badtimes,1)=swhistbins(LOCS(1));
 
 
%SWS time points
NREMtimes = (broadbandSlowWave >swthresh);


%% Then Divide EMG
numpeaks = 1;
numbins = 12;
if sum(isnan(EMG))>0
   error('EMG seems to have NaN values...') 
end

while numpeaks ~=2
    [EMGhist,EMGhistbins]= hist(EMG(NREMtimes==0),numbins);
    %[EMGhist,EMGhistbins]= hist(EMG,numbins);

    [PKS,LOCS] = findpeaks_SleepScore([0 EMGhist],'NPeaks',2);
    LOCS = sort(LOCS)-1;
    numbins = numbins+1;
    numpeaks = length(LOCS);
    
    if numpeaks ==100
        display('Something is wrong with your EMG')
        return
    end
end

betweenpeaks = EMGhistbins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks_SleepScore(-EMGhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

EMGthresh = betweenpeaks(diploc);

MOVtimes = (broadbandSlowWave<swthresh & EMG>EMGthresh);


%% Then Divide Theta
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

numbins = 12;
%numbins = 15; %for Poster...
while numpeaks ~=2 && numbins <=25
    [THhist,THhistbins]= hist(thratio(NREMtimes==0 & MOVtimes==0),numbins);

    [PKS,LOCS] = findpeaks_SleepScore(THhist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

if length(PKS)==2
    betweenpeaks = THhistbins(LOCS(1):LOCS(2));
    [dip,diploc] = findpeaks_SleepScore(-THhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

    THthresh = betweenpeaks(diploc);

    REMtimes = (broadbandSlowWave<swthresh & EMG<EMGthresh & thratio>THthresh);
else
    THthresh = 0;
%     REMtimes =(broadbandSlowWave<swthresh & EMG<EMGthresh);
end

histsandthreshs = v2struct(swhist,swhistbins,swthresh,EMGhist,EMGhistbins,EMGthresh,THhist,THhistbins,THthresh);

%% Ouput Structure: StateScoreMetrics
LFPparams = SleepScoreLFP.params;
THchanID = SleepScoreLFP.THchanID; SWchanID = SleepScoreLFP.SWchanID;

SleepScoreMetrics = v2struct(broadbandSlowWave,thratio,EMG,t_EMG,...
    t_clus,badtimes,reclength,histsandthreshs,LFPparams,THchanID,SWchanID,recordingname);
%save(matfilename,'SleepScoreMetrics');

StatePlotMaterials = v2struct(swFFTfreqs,swFFTspec,thFFTfreqs,thFFTspec);
%save(plotmaterialsfilename,'StatePlotMaterials'); 
