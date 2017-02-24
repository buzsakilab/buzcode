function [SWchannum,THchannum,swLFP,thLFP,t_LFP,Fs_save,SWfreqlist,SWweights] = PickSWTHChannel(datasetfolder,recordingname,figfolder,scoretime,SWWeightsName,Notch60Hz,NotchUnder3Hz,NotchHVS,NotchTheta,SWChannels,ThetaChannels);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%TO DO
%   -Change from GetLFP to LoadBinary or readmulti_ss
%% DEV
%datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/DTData/';
%recordingname = 'DT2_rPPC_rCCG_362um_218um_20160209_160209_183610';
% datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/~updated/Recordings (1)/';
% recordingname = 'c3po_160202';
% figfolder = '/Users/dlevenstein/Code Library/SleepScoreDevelopment/StateScoreFigures/';

%recname = 'c3po_160202';
%datasetfolder = '/Users/dlevenstein/Dropbox/Share Folders/Recordings/';

% datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/GGData/';
% recname = 'Rat08-20130717';

if ~exist('SWWeightsName','var')
    SWWeightsName = 'SWweights.mat';
end

xmlfilename = [datasetfolder,'/',recordingname,'/',recordingname,'.xml'];
if exist (fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.lfp']);
elseif exist (fullfile(datasetfolder,recordingname,[recordingname,'.eeg']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.eeg']);
else 
    display('No .lfp file')
end

%% FMA
% 
% SetCurrentSession(xmlfilename);
% global DATA
%nChannels = DATA.nChannels;

Par = LoadPar_SleepScore(xmlfilename);
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
numfreqs = 100;
freqlist = logspace(0,2,numfreqs);
window = 10;
noverlap = 9;
window = window*Fs;
noverlap = noverlap*Fs;


%% Pick channels to use
spkgroupchannels = [SpkGrps.Channels];

%Add reject channels here...
rejectchannels = [];
if exist(fullfile(datasetfolder,recordingname,'bad_channels.txt'),'file')%bad channels is an ascii/text file where all lines below the last blank line are assumed to each have a single entry of a number of a bad channel (base 0)
    t = ReadBadChannels_ss(fullfile(datasetfolder,recordingname));
    t = t+1;%account for offset
    rejectchannels = cat(1,rejectchannels(:),t(:));
end
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
allLFP = LoadBinary_Down_ss(rawlfppath,'frequency',Fs,...
    'nchannels',nChannels,'channels',usechannels+1,'downsample',downsamplefactor,...
    'start',scoretime(1),'duration',diff(scoretime));
Fs = Fs./downsamplefactor;
%+1 to channel number here to account for 0-indexing vs 1-indexing

%% For each channel, calculate the PC1 and check it
swhists = zeros(numhistbins,numSWChannels);
% pc1coeff = zeros(numfreqs,numusedchannels);
dipSW = zeros(numSWChannels,1);

THhist = zeros(numhistbins,numThetaChannels);
THmeanspec = zeros(numfreqs,numThetaChannels);
peakTH = zeros(numThetaChannels,1);

%% Get info to allow to pick SW channel
for idx = 1:numSWChannels;
%channum = 1;
    if mod(idx,10) == 1
        display(['Characterizing SW candidate channel ',num2str(idx),' of ',num2str(numSWChannels)])
    end

    tchanidx = SWChannels(idx);
    chanidx = find(usechannels==tchanidx);

    %% Get spectrogram
    %Calcualte Z-scored Spectrogram
    [FFTspec,FFTfreqs,t_FFT] = spectrogram(allLFP(:,chanidx),window,noverlap,freqlist,Fs);
    FFTspec = abs(FFTspec);
    [zFFTspec,mu,sig] = zscore(log10(FFTspec)');
    % Remove transients before calculating SW histogram
    %this should be it's own whole section - removing/detecting transients
    totz = zscore(abs(sum(zFFTspec')));
    badtimes = find(totz>5);
    zFFTspec(badtimes,:) = 0;
    
    %% PCA for Broadband Slow Wave
%     [COEFF, SCORE, LATENT] = pca(zFFTspec);
%    % broadbandSlowWave = SCORE(:,1);
%     
	%% Set Broadband filter weights for Slow Wave
    load(SWWeightsName)% 'SWweights.mat' by default
    assert(isequal(freqlist,SWfreqlist), 'spectrogram freqs.  are not what they should be...')
    
    %% Alter the filter weights if requested by the user
    if Notch60Hz
        SWweights(SWfreqlist<=62.5 & SWfreqlist>=57.5) = 0;
    end
    if NotchUnder3Hz
        SWweights(SWfreqlist<=3) = 0;
    end
    if NotchHVS
        SWweights(SWfreqlist<=18 & SWfreqlist>=12) = 0;
        SWweights(SWfreqlist<=10 & SWfreqlist>=4) = 0;
    end
    if NotchTheta
        SWweights(SWfreqlist<=10 & SWfreqlist>=4) = 0;
    end
    
    %% Calculate per-bin projections
    broadbandSlowWave = zFFTspec*SWweights';
    
    %% Smooth and 0-1 normalize
    smoothfact = 10; %units of si_FFT
    thsmoothfact = 10;
     
    broadbandSlowWave = smooth(broadbandSlowWave,smoothfact);
    broadbandSlowWave = (broadbandSlowWave-min(broadbandSlowWave))./max(broadbandSlowWave-min(broadbandSlowWave));

    %% Histogram and diptest of PC1
    histbins = linspace(0,1,numhistbins);
    [swhist]= hist(broadbandSlowWave,histbins);

    swhists(:,idx) = swhist;
%     pc1coeff(:,chanidx) = COEFF(:,1);
    dipSW(idx) = hartigansdiptest_ss(sort(broadbandSlowWave));
    
end

%% Get info to allow to pick Theta channel
for idx = 1:numThetaChannels;
%channum = 1;
    if mod(idx,10) == 1
        display(['Characterizing theta candidate channel ',num2str(idx),' of ',num2str(numSWChannels)])
    end

    tchanidx = ThetaChannels(idx);
    chanidx = find(usechannels==tchanidx);

    %% Get spectrogram
    %NarrowbandTheta
    %f_all = [3 16];
    f_all = [2 20];
    f_theta = [5 10];
    thfreqlist = logspace(log10(f_all(1)),log10(f_all(2)),numfreqs);

    [thFFTspec,thFFTfreqs] = spectrogram(allLFP(:,chanidx),window,noverlap,thfreqlist,Fs);
    thFFTspec = (abs(thFFTspec));

    thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
    thpower = sum((thFFTspec(thfreqs,:)),1);
    allpower = sum((thFFTspec),1);

    thratio = thpower./allpower;    %Narrowband Theta
    thratio = smooth(thratio,thsmoothfact);
    thratio = (thratio-min(thratio))./max(thratio-min(thratio));
    
    %% Histogram and diptest of Theta
    THhist(:,chanidx) = hist(thratio,histbins);
    dipTH(chanidx) = hartigansdiptest_ss(sort(thratio));
    
    %% Theta Peak in mean spectrum
    THmeanspec(:,idx) = (mean(thFFTspec,2));
    %THmeanspec(:,chanidx) = THmeanspec(:,chanidx)-min(THmeanspec(:,chanidx));
    meanthratio = sum((THmeanspec(thfreqs,idx)))./sum((THmeanspec(:,idx)));
    peakTH(idx) = meanthratio;
end

%% Sort by dip and pick channels
[~,dipsortSW] = sort(dipSW);
[~,dipsortTH] = sort(peakTH);

goodSWidx = dipsortSW(end);
goodTHidx = dipsortTH(end);

SWchannum = SWChannels(goodSWidx);
THchannum = ThetaChannels(goodTHidx);

downsample_save = 5;  %Not checked for bugs after adding...
Fs_save = Par.lfpSampleRate./downsample_save;
swthLFP = LoadBinary_Down_ss(rawlfppath,'frequency',Fs,...
    'downsample',downsample_save,...
    'nchannels',nChannels,'channels',[SWchannum+1,THchannum+1],...
    'start',scoretime(1),'duration',diff(scoretime));

swLFP = swthLFP(:,1);
thLFP = swthLFP(:,2);
t_LFP = [1:length(swLFP)]./Fs_save;

%% Find Inverted PC1s and flip them for plot
% invpc1 = mean(pc1coeff(freqlist<4,:))<0 & mean(pc1coeff(freqlist>50,:))>0;
% pc1coeff(:,invpc1) = -pc1coeff(:,invpc1);
% swhists(:,invpc1) = flipud(swhists(:,invpc1));
%% Test
%PC1 coefficients for NREM match
%Theta spectrum for isolated peak?


%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PC1 Weights and Coefficients


swfig = figure;
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

thfig = figure;
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
    [FFTspec,FFTfreqs,t_FFT] = spectrogram(allLFP(:,goodSWidx),window,noverlap,freqlist,Fs);
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

chanfig =figure;
	subplot(5,1,1:2)
        imagesc(t_FFT,log2(FFTfreqs),log10(FFTspec))
        axis xy
        LogScale_ss('y',2)
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        ylim([log2(FFTfreqs(1)) log2(FFTfreqs(end))+0.2])
        xlim(t_FFT([1,end]))
        ylabel({'LFP - FFT','f (Hz)'})
        title('SW Channel');
        
    subplot(5,1,3)
        plot(t_FFT,broadbandSlowWave,'k')
        xlim(t_FFT([1,end]))
        set(gca,'XTick',[]);
     
    %Calculate Theta ratio for plot/return    
    [thFFTspec,thFFTfreqs,t_FFT] = spectrogram(allLFP(:,goodTHidx),window,noverlap,thfreqlist,Fs);
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
        title('Theta Channel');
        xlim(t_FFT([1,end]))
        set(gca,'XTick',[]);
        
subplot(5,1,5)
        plot(t_FFT,thratio,'k')
        xlim(t_FFT([1,end]))
        set(gca,'XTick',[]);
        
saveas(chanfig,[figfolder,recordingname,'_SWTHChannels'],'jpeg')
%saveas(chanfig,[figfolder,recordingname,'_SWTHChannels'],'fig')
end

