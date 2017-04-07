function [broadbandSlowWave,thratio,EMG,t_EMG,t_FFT,badtimes, reclength,...
    histsandthreshs,FFTfreqs,FFTspec,thFFTfreqs,thFFTspec] = ClusterStates_GetParams(...
    swLFP,thLFP,EMG,sf_LFP,sf_EMG,figloc,recordingname,MinWinParams)
%StateID(LFP,thLFP,EMG,sf_LFP,sf_EMG,figloc,WSEpisodes)
%   Detailed explanation goes here
%
%
%
%Dependencies: IDXtoINT, INTtoIDX
%
%Last Updated: 1/31/16
%DLevenstein


%% Downsample and filter
%Make Downsample to niquest frequency

if sf_LFP == 1250
    downsamplefactor = 5;
elseif sf_LFP == 250
    downsamplefactor = 1;
elseif sf_LFP == 1000
    downsamplefactor = 4;
else
    display('sf not recognized... if only you made this able to set its own downsample...')
end
swLFP = downsample(swLFP,downsamplefactor);
thLFP = downsample(thLFP,downsamplefactor);
sf_LFP = sf_LFP/downsamplefactor;


%filtbounds = [0.5 120];
%display(['Filtering ',num2str(filtbounds(1)),'-',num2str(filtbounds(2)),' Hz...']);
%LFP = FiltNPhase(LFP, filtbounds, sf_LFP );


%% Calculate Spectrogram
%display('FFT Spectrum for Broadband LFP')

freqlist = logspace(0,2,100);
%freqlist = linspace(0.5,55.5,111);
window = 10;   %s
noverlap = 9;  %s
window = window*sf_LFP;
noverlap = noverlap*sf_LFP;
[FFTspec,FFTfreqs,t_FFT] = spectrogram(swLFP,window,noverlap,freqlist,sf_LFP);
FFTspec = abs(FFTspec);
[zFFTspec,mu,sig] = zscore(log10(FFTspec)');

    %% Remove transients before calculating SW histogram
    %this should be it's own whole section - removing/detecting transients
totz = zscore(abs(sum(zFFTspec')));
badtimes = find(totz>5);
zFFTspec(badtimes,:) = 0;

%% PCA for Broadband Slow Wave
%  [COEFF, SCORE, ~, ~, EXPLAINED] = pca(zFFTspec);
%  % broadbandSlowWave = SCORE(:,1);
% PC1weights = COEFF(:,1);
% PC1expvar = EXPLAINED(1);
 
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


[thFFTspec,thFFTfreqs] = spectrogram(thLFP,window,noverlap,freqlist,sf_LFP);
thFFTspec = (abs(thFFTspec));
[~,mu_th,sig_th] = zscore(log10(thFFTspec)');

thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
allpower = sum((thFFTspec),1);
thpower = sum((thFFTspec(thfreqs,:)),1);

thratio = thpower./allpower;    %Narrowband Theta
thratio = smooth(thratio,thsmoothfact);
thratio = (thratio-min(thratio))./max(thratio-min(thratio));
 
%% EMG
dtEMG = 1/sf_EMG;
t_EMG = (1:length(EMG))*dtEMG;
EMG = smooth(EMG,smoothfact/dtEMG);
EMG = (EMG-min(EMG))./max(EMG-min(EMG));

reclength = round(t_EMG(end));

%downsample to FFT time points;
[~,t_intersect] = intersect(t_EMG,t_FFT);
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


