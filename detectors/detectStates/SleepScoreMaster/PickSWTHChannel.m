function [SleepScoreLFP] = PickSWTHChannel(basePath,figfolder,scoretime,SWWeightsName,Notch60Hz,NotchUnder3Hz,NotchHVS,NotchTheta,SWChannels,ThetaChannels,rejectchannels,OVERWRITE);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%% Buzcode name of the SleepScoreLFP.LFP.mat file
[datasetfolder,recordingname,extension] = fileparts(basePath);
recordingname = [recordingname extension];

matfilename = fullfile(basePath,[recordingname,'.SleepScoreLFP.LFP.mat']);

saveFiles = true;
%% Check if SleepScoreLFP has already been claculated for this recording
%If the SleepScoreLFP file already exists, load and return with SleepScoreLFP in hand
if exist(matfilename,'file') && ~OVERWRITE
    display('SleepScoreLFP already calculated - loading from SleepScoreLFP.LFP.mat')
    load(matfilename)
    if ~exist('SleepScoreLFP','var')
        display([matfilename,' does not contain a variable called SleepScoreLFP'])
    end
    return
end
display('Picking SW and TH Channels for SleepScoreLFP.LFP.mat')

%%

xmlfilename = [datasetfolder,'/',recordingname,'/',recordingname,'.xml'];
if exist (fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.lfp']);
elseif exist (fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.lfp']);
elseif exist (fullfile(datasetfolder,recordingname,[recordingname,'.eeg']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.eeg']);
else 
    display('No .lfp file')
end

%% FMA
Par = bz_getSessionInfo(basePath);
Fs = Par.lfpSampleRate; % Hz, LFP sampling rate
nChannels = Par.nChannels;

if isfield(Par,'SpkGrps')
    SpkGrps = Par.SpkGrps;
elseif isfield(Par,'AnatGrps')
    SpkGrps = Par.AnatGrps;
    display('No SpikeGroups, Using AnatomyGroups')
else
    display('No SpikeGroups...')
end

%% Hist/Freqs Parms
numhistbins = 21;
histbins = linspace(0,1,numhistbins);
numfreqs = 100;
swFFTfreqs = logspace(0,2,numfreqs);
window = 10;
noverlap = 9;
window = window*Fs;
noverlap = noverlap*Fs;

%Smoothing Parameters
smoothfact = 10; %units of si_FFT
thsmoothfact = 10;

%For SW calculation
%Load the slowwave filter weights
if ~exist('SWWeightsName','var')
    SWWeightsName = 'SWweights.mat';
end
load(SWWeightsName)% 'SWweights.mat' by default
%Alter the filter weights if requested by the user
if Notch60Hz; SWweights(SWfreqlist<=62.5 & SWfreqlist>=57.5) = 0; end
if NotchUnder3Hz; SWweights(SWfreqlist<=3) = 0; end
if NotchHVS
    SWweights(SWfreqlist<=18 & SWfreqlist>=12) = 0;
    SWweights(SWfreqlist<=10 & SWfreqlist>=4) = 0;
end
if NotchTheta; SWweights(SWfreqlist<=10 & SWfreqlist>=4) = 0; end

assert(isequal(swFFTfreqs,SWfreqlist), 'spectrogram freqs.  are not what they should be...')
   

%For Theta Calculation
f_all = [2 20];
f_theta = [5 10];
thFFTfreqs = logspace(log10(f_all(1)),log10(f_all(2)),numfreqs);


%% Pick channels to use
spkgroupchannels = [SpkGrps.Channels];

if sum(SWChannels)>0 && sum(ThetaChannels)>0%use all channels unless SWChannels and ThetaChannels are specified... if both specified then we know those are the only good ones
    goodchannels = union(SWChannels,ThetaChannels);
    badchannels = setdiff(spkgroupchannels,goodchannels);
    rejectchannels = union(rejectchannels,badchannels);
end

usechannels = setdiff(spkgroupchannels,rejectchannels);
numusedchannels = length(usechannels);

%% Handle specific candidacy of certain channels for SW vs Theta
if sum(SWChannels)==0
    SWChannels = usechannels;
else
    SWChannels = setdiff(SWChannels,rejectchannels);
end
numSWChannels = length(SWChannels);

if sum(ThetaChannels)==0
    ThetaChannels = usechannels;
else
    ThetaChannels = setdiff(ThetaChannels,rejectchannels);
end
numThetaChannels = length(ThetaChannels);

%% Load LFP files from .lfp
downsamplefactor = 10;
allLFP = bz_LoadBinary(rawlfppath,'frequency',Fs,...
    'nchannels',nChannels,'channels',usechannels+1,'downsample',downsamplefactor,...
    'start',scoretime(1),'duration',diff(scoretime));
%allLFP = double(allLFP); % hack fix
Fs = Fs./downsamplefactor;

%% For each channel, calculate the PC1 and check it
swhists = zeros(numhistbins,numSWChannels);
dipSW = zeros(numSWChannels,1);

THhist = zeros(numhistbins,numThetaChannels);
THmeanspec = zeros(numfreqs,numThetaChannels);
peakTH = zeros(numThetaChannels,1);


%% Get info to allow to pick SW channel
parfor_progress(numSWChannels);
tstart = tic;
parfor idx = 1:numSWChannels;
%channum = 1;

    %Progress Counter
    timespent=toc(tstart);
    percdone = parfor_progress;
    
    estimatedtotal = timespent./(percdone./100);
    estimatedremaining = estimatedtotal-timespent;
   %if mod(idx,10) == 1
   %fprintf('\r'); % delete previous counter display
        display(['SW Channels - Percent Complete: ',num2str(round(percdone)),...
            '.  Time Spent: ',num2str(round(timespent./60)),...
            '.  Est. Total Time: ',num2str(round(estimatedtotal./60)),...
            'min.  ETR: ',num2str(round(estimatedremaining./60)),'min.'])
  % end

    %% Get spectrogram
    %Calcualte Z-scored Spectrogram
    LFPchanidx = find(usechannels==SWChannels(idx));
    FFTspec = spectrogram(single(allLFP(:,LFPchanidx)),window,noverlap,swFFTfreqs,Fs);
    FFTspec = abs(FFTspec);
    [zFFTspec,mu,sig] = zscore(log10(FFTspec)');
    % Remove transients before calculating SW histogram
    %this should be it's own whole section - removing/detecting transients
    totz = zscore(abs(sum(zFFTspec')));
    badtimes = find(totz>5);
    zFFTspec(badtimes,:) = 0;
  
    %% Calculate per-bin weights onto SlowWave
    broadbandSlowWave = zFFTspec*SWweights';
    broadbandSlowWave = smooth(broadbandSlowWave,smoothfact);
    broadbandSlowWave = (broadbandSlowWave-min(broadbandSlowWave))./max(broadbandSlowWave-min(broadbandSlowWave));

    %% Histogram and diptest of Slow Wave Power
    [swhist]= hist(broadbandSlowWave,histbins);
    
    %Record the histogram and dip score for later comparison between chans
    swhists(:,idx) = swhist;
    dipSW(idx) = hartigansdiptest_ss(sort(broadbandSlowWave));
end
parfor_progress(0);
%% Get info to allow to pick Theta channel
parfor_progress(numSWChannels);
tstart = tic;
parfor idx = 1:numThetaChannels;
%channum = 1;
    %Progress Counter
    timespent=toc(tstart);
    percdone = parfor_progress;
    
    estimatedtotal = timespent./(percdone./100);
    estimatedremaining = estimatedtotal-timespent;
   %if mod(idx,10) == 1
   %fprintf('\r'); % delete previous counter display
        display(['TH Channels - Percent Complete: ',num2str(round(percdone)),...
            '.  Time Spent: ',num2str(round(timespent./60)),...
            '.  Est. Total Time: ',num2str(round(estimatedtotal./60)),...
            'min.  ETR: ',num2str(round(estimatedremaining./60)),'min.'])
  % end

    %% Get spectrogram and calculate theta ratio
    LFPchanidx = find(usechannels==ThetaChannels(idx));
    thFFTspec = spectrogram(single(allLFP(:,LFPchanidx)),window,noverlap,thFFTfreqs,Fs);
    thFFTspec = (abs(thFFTspec));

    thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
    thpower = sum((thFFTspec(thfreqs,:)),1);
    allpower = sum((thFFTspec),1);

    thratio = thpower./allpower;    %Narrowband Theta
    thratio = smooth(thratio,thsmoothfact);
    thratio = (thratio-min(thratio))./max(thratio-min(thratio));
    
    %% Histogram and diptest of Theta
    THhist(:,idx) = hist(thratio,histbins);
    %Dip test of theta doesn't get used... could be incorporated for
    %selection?
    %dipTH(idx) = hartigansdiptest_ss(sort(thratio));
    
    %% Ratio of Theta Peak to sorrounding in mean spectrum (for selection)
    meanspec = (mean(thFFTspec,2));
    meanthratio = sum((meanspec(thfreqs)))./sum((meanspec(:)));
    
    %Record the spec and peak ratio for later comparison between chans
    THmeanspec(:,idx) = meanspec;
    peakTH(idx) = meanthratio;
end
parfor_progress(0);
%% Sort by dip (bimodality) and pick channels
[~,dipsortSW] = sort(dipSW);
[~,dipsortTH] = sort(peakTH);

goodSWidx = dipsortSW(end); %Channel list Index of the 
goodTHidx = dipsortTH(end); %best SW and theta channels

SWchanID = SWChannels(goodSWidx);      %Channel IDnumber of the 
THchanID = ThetaChannels(goodTHidx);   %best SW and theta channels

%% Load the best channels at sampling frequency needed for clustering later
downsample_save = Par.lfpSampleRate./250;
sf = Par.lfpSampleRate./downsample_save;
swthLFP = bz_LoadBinary(rawlfppath,'frequency',Par.lfpSampleRate,...
    'downsample',downsample_save,...
    'nchannels',nChannels,'channels',[SWchanID,THchanID]+1,...
    'start',scoretime(1),'duration',diff(scoretime));

swLFP = (swthLFP(:,1));
thLFP = (swthLFP(:,2));
t = [1:length(swLFP)]./sf;


%% SleepScoreLFP output

params = v2struct(SWfreqlist,SWweights,SWWeightsName,Notch60Hz,...
    NotchUnder3Hz,NotchHVS,NotchTheta);
    
SleepScoreLFP = v2struct(thLFP,swLFP,THchanID,SWchanID,sf,t,params);


if saveFiles
    %Need to update to Save in buzcode format for lfp.mat
    save(matfilename,'SleepScoreLFP');
end


%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PC1 Weights and Coefficients


swfig = figure('visible','off');
%     subplot(2,2,1)
%         imagesc(log2(FFTfreqs),1:numusedchannels,pc1coeff(:,dipsortSW)')
%         ylabel('Channel #');xlabel('f (Hz)')
%         LogScale_ss('x',2)
%         axis xy
%         title('PC1 Frequency Coefficients: All Channels') 
    subplot(2,2,2)
        imagesc(histbins,1:numusedchannels,swhists(:,dipsortSW)')
        ylabel('Channel #');xlabel('SW Band projection weight')
        title('SW Band Projection Histogram: All Channels')
        axis xy
%     subplot(2,2,3)
%         set(gca,'ColorOrder',RainbowColors_ss(length(dipsortSW)))
%         hold all
%         plot(log2(FFTfreqs),pc1coeff')  
%         plot(log2(FFTfreqs),pc1coeff(:,goodSWidx)','k','LineWidth',1)
%         plot(log2(FFTfreqs([1 end])),[0 0],'k')
%         ylabel('PC1 Coefficient');xlabel('f (Hz)')
%         LogScale_ss('x',2)
%         title('PC1 Frequency Coefficients: All Channels')
    subplot(2,2,4)
        set(gca,'ColorOrder',RainbowColors_ss(length(dipsortSW)))
        hold all
        plot(histbins,swhists')
        plot(histbins,swhists(:,goodSWidx)','k','LineWidth',1)
        ylabel('hist');xlabel('SW Band projection weight')
        title('SW Band Projection Histogram: All Channels')
        
saveas(swfig,[figfolder,recordingname,'_FindBestSW'],'jpeg')

%% Theta Hist and Coefficients

thfig = figure('visible','off');
    subplot(2,2,1)
        imagesc(log2(thFFTfreqs),1:numusedchannels,THmeanspec(:,dipsortTH)')
        ylabel('Channel #');xlabel('f (Hz)')
        LogScale_ss('x',2)
        axis xy
        title('Spectrum: All Channels') 
    subplot(2,2,2)
        imagesc(histbins,1:numusedchannels,THhist(:,dipsortTH)')
        ylabel('Channel #');xlabel('Theta projection weight')
        title('Theta Ratio Histogram: All Channels')
        axis xy
    subplot(2,2,3)
        set(gca,'ColorOrder',RainbowColors_ss(length(dipsortTH)))
        hold all
        plot(log2(thFFTfreqs),THmeanspec')  
        plot(log2(thFFTfreqs),THmeanspec(:,goodTHidx)','k','LineWidth',1)
        plot(log2(f_theta(1))*[1 1],get(gca,'ylim'),'k')
        plot(log2(f_theta(2))*[1 1],get(gca,'ylim'),'k')
        ylabel('Power');xlabel('f (Hz)')
        xlim(log2(thFFTfreqs([1 end])))
        LogScale_ss('x',2)
        title('Spectrum: All Channels')
    subplot(2,2,4)
        set(gca,'ColorOrder',RainbowColors_ss(length(dipsortTH)))
        hold all
        plot(histbins,THhist')
        plot(histbins,THhist(:,goodTHidx)','k','LineWidth',1)
        ylabel('hist');xlabel('Theta projection weight')
        title('Theta Ratio Histogram: All Channels') 
        
saveas(thfig,[figfolder,recordingname,'_FindBestTH'],'jpeg')
%% Show Channels


    %Calculate PC1 for plot/return
    [FFTspec,swFFTfreqs,t_FFT] = spectrogram(double(allLFP(:,goodSWidx)),window,noverlap,swFFTfreqs,Fs);
    FFTspec = abs(FFTspec);
    [zFFTspec,mu,sig] = zscore(log10(FFTspec)');

    totz = zscore(abs(sum(zFFTspec')));
    badtimes = find(totz>5);
    zFFTspec(badtimes,:) = 0;
    
     %[COEFF, SCORE, LATENT] = pca(zFFTspec);
    %broadbandSlowWave = SCORE(:,1);
     broadbandSlowWave = zFFTspec*SWweights';
     broadbandSlowWave = smooth(broadbandSlowWave,smoothfact);
     broadbandSlowWave = (broadbandSlowWave-min(broadbandSlowWave))./max(broadbandSlowWave-min(broadbandSlowWave));

chanfig =figure('visible','off');
	subplot(5,1,1:2)
        imagesc(t_FFT,log2(swFFTfreqs),log10(FFTspec))
        axis xy
        LogScale_ss('y',2)
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        ylim([log2(swFFTfreqs(1)) log2(swFFTfreqs(end))+0.2])
        xlim(t_FFT([1,end]))
        ylabel({'LFP - FFT','f (Hz)'})
        title(['SW Channel:',num2str(SWchanID)]);
        
    subplot(5,1,3)
        plot(t_FFT,broadbandSlowWave,'k')
        xlim(t_FFT([1,end]))
        set(gca,'XTick',[]);
     
    %Calculate Theta ratio for plot/return    
    [thFFTspec,thFFTfreqs,t_FFT] = spectrogram(double(allLFP(:,goodTHidx)),window,noverlap,thFFTfreqs,Fs);
    thFFTspec = abs(thFFTspec);
    [zFFTspec,mu,sig] = zscore(log10(thFFTspec)');
        
    thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
    thpower = sum((thFFTspec(thfreqs,:)),1);
    allpower = sum((thFFTspec),1);

    thratio = thpower./allpower;    %Narrowband Theta
    thratio = smooth(thratio,thsmoothfact);
    thratio = (thratio-min(thratio))./max(thratio-min(thratio));
    
subplot(5,1,4)
 %   plot(allLFP(:,1),allLFP(:,goodSWidx),'k')
    
        imagesc(t_FFT,log2(thFFTfreqs),log10(thFFTspec))
        hold on
        plot(t_FFT([1,end]),log2(f_theta([1,1])),'w')
        plot(t_FFT([1,end]),log2(f_theta([2,2])),'w')
        axis xy
        LogScale_ss('y',2)
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        ylim([log2(thFFTfreqs(1)) log2(thFFTfreqs(end))+0.2])
        ylabel({'LFP - FFT','f (Hz)'})
        title(['Theta Channel: ',num2str(THchanID)]);
        xlim(t_FFT([1,end]))
        set(gca,'XTick',[]);
        
subplot(5,1,5)
        plot(t_FFT,thratio,'k')
        xlim(t_FFT([1,end]))
        set(gca,'XTick',[]);
        
saveas(chanfig,[figfolder,recordingname,'_SWTHChannels'],'jpeg')
%saveas(chanfig,[figfolder,recordingname,'_SWTHChannels'],'fig')
end

