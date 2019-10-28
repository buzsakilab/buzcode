function [SleepScoreLFP,PickChannelStats] = PickSWTHChannel(basePath,scoretime,SWWeightsName,...
    Notch60Hz,NotchUnder3Hz,NotchHVS,NotchTheta,SWChannels,ThetaChannels,...
    rejectchannels,OVERWRITE,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Optional
%       'noPrompts'     (default: false) prevents prompts about saving/adding metadata
%       'ignoretime'    time windows to ignore when claculting SW/TH for
%                       channel selection
%
%% Buzcode name of the SleepScoreLFP.LFP.mat file
p = inputParser;
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'saveFiles',false,@islogical);
addParameter(p,'ignoretime',[]);
addParameter(p,'window',2);
addParameter(p,'SHOWFIG',false);
addParameter(p,'smoothfact',15);
addParameter(p,'IRASA',true);
addParameter(p,'downsamplefactor',5);

parse(p,varargin{:})
noPrompts = p.Results.noPrompts;
saveFiles = p.Results.saveFiles;
ignoretime = p.Results.ignoretime;
window = p.Results.window; 
smoothfact = p.Results.smoothfact; 
IRASA = p.Results.IRASA; 
SHOWFIG = p.Results.SHOWFIG; 
downsamplefactor = p.Results.downsamplefactor; 


[datasetfolder,recordingname,extension] = fileparts(basePath);
recordingname = [recordingname extension];

matfilename = fullfile(basePath,[recordingname,'.SleepScoreLFP.LFP.mat']);

figfolder = [fullfile(basePath,'StateScoreFigures'),'/'];
if ~exist(figfolder,'dir') & saveFiles
    mkdir(figfolder)
end

%% Hist/Freqs Parms
numhistbins = 21;
histbins = linspace(0,1,numhistbins);
numfreqs = 100;
swFFTfreqs = logspace(0,2,numfreqs);
noverlap = 0; %Updated to speed up (don't need to sample at fine time resolution for channel selection)


%For SW calculation
%Load the slowwave filter weights
if ~exist('SWWeightsName','var') | strcmp(SWWeightsName,'PSS')
    %SWWeightsName = 'SWweights.mat';
    SWweights = 'PSS';
    SWWeightsName = 'PSS';
    SWfreqlist = logspace(0.5,2,numfreqs); %should get this from bz_PowerSpectrumSlope...
else
    load(SWWeightsName)% 'SWweights.mat' by default
    %Alter the filter weights if requested by the user
    if Notch60Hz; SWweights(SWfreqlist<=62.5 & SWfreqlist>=57.5) = 0; end
    if NotchUnder3Hz; SWweights(SWfreqlist<=3) = 0; end
    if NotchHVS
        SWweights(SWfreqlist<=18 & SWfreqlist>=12) = 0;
        SWweights(SWfreqlist<=10 & SWfreqlist>=4) = 0;
    end
    if NotchTheta; SWweights(SWfreqlist<=10 & SWfreqlist>=4) = 0; end

    assert(isequal(swFFTfreqs,SWfreqlist),...
        'spectrogram freqs.  are not what they should be...')
end 

%For Theta Calculation
f_all = [2 20];
f_theta = [5 10];
thFFTfreqs = logspace(log10(f_all(1)),log10(f_all(2)),numfreqs);

%% Check if SleepScoreLFP has already been claculated for this recording
%If the SleepScoreLFP file already exists, load and return with SleepScoreLFP in hand
if exist(matfilename,'file') && ~OVERWRITE
    display('SleepScoreLFP already calculated - loading from SleepScoreLFP.LFP.mat')
    load(matfilename)
    if ~exist('SleepScoreLFP','var')
        display([matfilename,' does not contain a variable called SleepScoreLFP'])
    end
    
    if ~isequal(SleepScoreLFP.params.SWWeightsName,SWWeightsName)
       display(['SlowWave Method used for Channel selection doesn''t match, updating to ',SWWeightsName])
       SleepScoreLFP.params.SWfreqlist_selection = SleepScoreLFP.params.SWfreqlist;
       SleepScoreLFP.params.SWweights_selection = SleepScoreLFP.params.SWweights;
       SleepScoreLFP.params.SWWeightsName_selection = SleepScoreLFP.params.SWWeightsName;
       SleepScoreLFP.params.SWfreqlist = SWfreqlist;
       SleepScoreLFP.params.SWweights = SWweights;
       SleepScoreLFP.params.SWWeightsName = SWWeightsName;
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

%% Pick channels to use

Par = bz_getSessionInfo(basePath,'noPrompts',noPrompts);
nChannels = Par.nChannels;

%Remove spike groups requirement DL9/3/19
% if isfield(Par,'SpkGrps')
%     SpkGrps = Par.SpkGrps;
% elseif isfield(Par,'AnatGrps')
%     SpkGrps = Par.AnatGrps;
%     display('No SpikeGroups, Using AnatomyGroups')
% else
%     display('No SpikeGroups...')
% end

% spkgroupchannels = [SpkGrps.Channels];
% 
% try %In case some channels are in AnatGrps but not SpkGrps
% anatgoupchannels = [Par.AnatGrps.Channels];
% spkgroupchannels = union(spkgroupchannels,anatgoupchannels);
% catch
% end

if sum(SWChannels)>0 && sum(ThetaChannels)>0%use all channels unless SWChannels and ThetaChannels are specified... if both specified then we know those are the only good ones
    goodchannels = union(SWChannels,ThetaChannels);
    badchannels = setdiff(Par.channels,goodchannels);
    rejectchannels = union(rejectchannels,badchannels);
end

usechannels = setdiff(Par.channels,rejectchannels);
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
allLFP = bz_GetLFP(usechannels,'basepath',basePath,'basename',recordingname,...
    'downsample',downsamplefactor,'intervals',scoretime,'noPrompts',noPrompts);
Fs = allLFP.samplingRate;


%% Set up containers for parfor loop
swhists = zeros(numhistbins,numSWChannels);
dipSW = zeros(numSWChannels,1);

THhist = zeros(numhistbins,numThetaChannels);
THmeanspec = zeros(numfreqs,numThetaChannels);
peakTH = zeros(numThetaChannels,1);

%% Get info to allow to pick SW channel
%parfor_progress(numSWChannels);
%tstart = tic;
parfor idx = 1:numSWChannels;
    %Progress Counter
%     timespent=toc(tstart);
%     percdone = parfor_progress;
%     
%     estimatedtotal = timespent./(percdone./100);
%     estimatedremaining = estimatedtotal-timespent;
   %if mod(idx,10) == 1
   %fprintf('\r'); % delete previous counter display
%         display(['SW Channels - Percent Complete: ',num2str(round(percdone)),...
%             '.  Time Spent: ',num2str(round(timespent./60)),...
%             '.  Est. Total Time: ',num2str(round(estimatedtotal./60)),...
%             'min.  ETR: ',num2str(round(estimatedremaining./60)),'min.'])
  % end
    bz_Counter(idx,numSWChannels,'SW Channels')
    %% Get Slow Wave signal from weighted or slope of the spectrogram
    %Calcualte Z-scored Spectrogram
    LFPchanidx = find(usechannels==SWChannels(idx));
    
    if strcmp(SWweights,'PSS')
        [specslope,~] = bz_PowerSpectrumSlope(allLFP,window,window-noverlap,...
            'channels',SWChannels(idx),'frange',[4 90],'nfreqs',100,'IRASA',IRASA);
        broadbandSlowWave = specslope.data;
        SWfreqlist(idx,:) = specslope.freqs;
        specdt = 1./specslope.samplingRate;
        t_FFT = specslope.timestamps;

    else
        [FFTspec,~,t_FFT] = spectrogram(single(allLFP.data(:,LFPchanidx)),window*Fs,noverlap*Fs,swFFTfreqs,Fs);
        t_FFT = t_FFT+allLFP.timestamps(1); %Offset for scoretime start
        FFTspec = abs(FFTspec);
        [zFFTspec,mu,sig] = zscore(log10(FFTspec)');
        % Remove transients before calculating SW histogram
        %this should be it's own whole section - removing/detecting transients
        totz = zscore(abs(sum(zFFTspec')));
        badtimes = find(totz>5);
        zFFTspec(badtimes,:) = 0;
        
        specdt = mode(diff(t_FFT));
        %Calculate per-bin weights onto SlowWave
        broadbandSlowWave = zFFTspec*SWweights';
        
    end
    
    
    broadbandSlowWave = smooth(broadbandSlowWave,smoothfact./specdt);
    
    %Remove ignoretimes (after smoothing), before normalizoing
    if ~isempty(ignoretime)
        ignoretimeIDX = InIntervals(t_FFT,ignoretime);
        broadbandSlowWave(ignoretimeIDX) = [];  %% this introduces gaps in a continuous signal... maybe nan fill?
    end
    broadbandSlowWave = bz_NormToRange(broadbandSlowWave,[0 1]);

    %% Histogram and diptest of Slow Wave Power
    [swhist]= hist(broadbandSlowWave,histbins);
    
    %Record the histogram and dip score for later comparison between chans
    swhists(:,idx) = swhist;
    dipSW(idx) = hartigansdiptest_ss(sort(broadbandSlowWave));
end
SWfreqlist = SWfreqlist(1,:);
%parfor_progress(0);

%% Get info to allow to pick Theta channel
%parfor_progress(numSWChannels);
%tstart = tic;
parfor idx = 1:numThetaChannels;
%channum = 1;
    %Progress Counter
    bz_Counter(idx,numThetaChannels,'TH Channels')

    if IRASA && strcmp(SWweights,'PSS') %(new way... peak above 1/f)
        
        [specslope,~] = bz_PowerSpectrumSlope(allLFP,window,window-noverlap,...
            'channels',ThetaChannels(idx),'frange',f_all,'nfreqs',100,'IRASA',IRASA);
        specdt = 1./specslope.samplingRate;
        t_FFT = specslope.timestamps;
        thFFTspec = specslope.resid';
        thFFTspec(thFFTspec<0)=0;
        
        % Remove transients before calculating SW histogram
        %zFFTspec = NormToInt(spec.amp,'modZ');
        %totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
        %badtimes_TH = find(totz>3);
        
        thFFTfreqs = specslope.freqs';
        thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
        thratio = max((thFFTspec(thfreqs,:)),[],1);
        thratio = smooth(thratio,smoothfact./specdt);

    else %(old way... ratio)
        %% Get spectrogram and calculate theta ratio
        %HERE: MATCH from GetMetrics for IRASA method
        LFPchanidx = find(usechannels==ThetaChannels(idx));
        [thFFTspec,~,t_FFT] = spectrogram(single(allLFP.data(:,LFPchanidx)),window*Fs,noverlap*Fs,thFFTfreqs,Fs);
        t_FFT = t_FFT+allLFP.timestamps(1); %Offset for scoretime start
        specdt = mode(diff(t_FFT));
        thFFTspec = (abs(thFFTspec));

        thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
        thpower = sum((thFFTspec(thfreqs,:)),1);
        allpower = sum((thFFTspec),1);

        thratio = thpower./allpower;    %Narrowband Theta
        thratio = smooth(thratio,smoothfact./specdt);


    end
    
    %Remove ignoretimes (after smoothing), before normalizoing
    if ~isempty(ignoretime)
        ignoretimeIDX = InIntervals(t_FFT,ignoretime);
        thratio(ignoretimeIDX) = [];
    end
    thratio = bz_NormToRange(thratio,[0 1]);
    %Histogram and diptest of Theta
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
%parfor_progress(0);
%% Sort by dip (bimodality) and pick channels
[~,dipsortSW] = sort(dipSW);
[~,dipsortTH] = sort(peakTH);

goodSWidx = dipsortSW(end); %Channel list Index of the 
goodTHidx = dipsortTH(end); %best SW and theta channels

SWchanID = SWChannels(goodSWidx);      %Channel IDnumber of the 
THchanID = ThetaChannels(goodTHidx);   %best SW and theta channels

%% Load the best channels at sampling frequency needed for clustering later
downsample_save = Par.lfpSampleRate./250;
swthLFP = bz_GetLFP([SWchanID,THchanID],'basepath',basePath,'basename',recordingname,...
    'downsample',downsample_save,'intervals',scoretime,'noPrompts',noPrompts);

swLFP = (swthLFP.data(:,1));
thLFP = (swthLFP.data(:,2));
t = swthLFP.timestamps;
sf = swthLFP.samplingRate;


%% SleepScoreLFP output

PickChannelStats = v2struct(dipSW,peakTH,SWChannels,ThetaChannels,...
    swhists,THmeanspec);

params = v2struct(SWfreqlist,SWweights,SWWeightsName,Notch60Hz,...
    NotchUnder3Hz,NotchHVS,NotchTheta,ignoretime,window,smoothfact,IRASA);
    
SleepScoreLFP = v2struct(thLFP,swLFP,THchanID,SWchanID,sf,t,params);

try
    bz_tagChannel(basePath,SWchanID,'NREMDetectionChan','noPrompts',true);
    bz_tagChannel(basePath,THchanID,'ThetaChan','noPrompts',true);
catch
    display('Unable to save channel tags in sessionInfo')
end

if saveFiles
    %Need to update to Save in buzcode format for lfp.mat
    save(matfilename,'SleepScoreLFP');
end


%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Old theta distribution panel
% 
%     subplot(2,2,2)
%         imagesc(histbins,1:numusedchannels,THhist(:,dipsortTH)')
%         ylabel('Channel #');xlabel('Theta projection weight')
%         title('Theta Ratio Histogram: All Channels')
%         axis xy
% 
%     subplot(2,2,4)
%         set(gca,'ColorOrder',RainbowColors_ss(length(dipsortTH)))
%         hold all
%         plot(histbins,THhist')
%         plot(histbins,THhist(:,goodTHidx)','k','LineWidth',1)
%         ylabel('hist');xlabel('Theta projection weight')
%         title('Theta Ratio Histogram: All Channels') 
        
%% Show Channels
if SHOWFIG
    chanfig =figure('visible','on');
else
    chanfig =figure('visible','off');
end
normhists = bz_NormToRange(swhists,[0 numusedchannels.*0.6]);
normTHspec = bz_NormToRange(THmeanspec,[0 numusedchannels.*0.6]);
    subplot(5,2,[7 9])
        imagesc(histbins,1:numusedchannels,swhists(:,dipsortSW)')
        axis xy
        hold on
        plot(histbins,normhists',...
            'color',0.9.*[1 1 1],'linewidth',0.1)
        plot(histbins,normhists(:,goodSWidx)','k','LineWidth',1)
        ylabel('Channel #');xlabel('SW weight')
        title('SW Histogram: All Channels')

    subplot(5,2,[8 10])
        imagesc(log2(thFFTfreqs),1:numusedchannels,THmeanspec(:,dipsortTH)')
        axis xy
        hold on
        plot(log2(thFFTfreqs),normTHspec',...
            'color',0.9.*[1 1 1],'linewidth',0.1)
        plot(log2(thFFTfreqs),normTHspec(:,goodTHidx)','k','LineWidth',1)
        ylabel('Channel #');xlabel('f (Hz)')
        LogScale_ss('x',2)
        axis xy
        title('TH Spectrum: All Channels') 

    %Calculate Slow Wave for figure
    if strcmp(SWweights,'PSS')
        [specslope,spec] = bz_PowerSpectrumSlope(allLFP,window,window-noverlap,...
            'channels',SWChannels(goodSWidx),'frange',[4 90]);
        broadbandSlowWave = specslope.data;
        specdt = 1./specslope.samplingRate;
        t_FFT = specslope.timestamps;
        t_spec = spec.timestamps;
        FFTspec = spec.amp';
        swFFTfreqs = spec.freqs;
        [zFFTspec,mu,sig] = zscore((FFTspec));

    else
        [FFTspec,~,t_FFT] = spectrogram(single(allLFP.data(:,goodSWidx)),window*Fs,noverlap*Fs,swFFTfreqs,Fs);
        t_FFT = t_FFT+allLFP.timestamps(1); %Offset for scoretime start
        t_spec = t_FFT;
        FFTspec = abs(FFTspec);
        [zFFTspec,mu,sig] = zscore(log10(FFTspec)');
        % Remove transients before calculating SW histogram
        %this should be it's own whole section - removing/detecting transients
        totz = zscore(abs(sum(zFFTspec')));
        badtimes = find(totz>5);
        zFFTspec(badtimes,:) = 0;
        
        specdt = mode(diff(t_FFT));
        %Calculate per-bin weights onto SlowWave
        broadbandSlowWave = zFFTspec*SWweights';
        
    end
 
    broadbandSlowWave = smooth(broadbandSlowWave,smoothfact./specdt);
    %Remove ignoretimes (after smoothing), before normalizoing
    if ~isempty(ignoretime)
        ignoretimeIDX = InIntervals(t_FFT,ignoretime);
        broadbandSlowWave(ignoretimeIDX) = [];
        t_FFT(ignoretimeIDX) = [];
    end
     
	subplot(5,1,1:2)
        imagesc(t_spec,log2(swFFTfreqs),(FFTspec))
        axis xy; hold on
        plot(t_FFT,bz_NormToRange(broadbandSlowWave,log2(swFFTfreqs([1 end]))),'k','Linewidth',0.1)
        LogScale_ss('y',2)
        caxis([min(mu)-2*max(sig) max(mu)+2*max(sig)])
        ylim([log2(swFFTfreqs(1)) log2(swFFTfreqs(end))+0.2])
        xlim(t_FFT([1,end]))
        
        ylabel({'LFP - FFT','f (Hz)'})
        title(['SW Channel:',num2str(SWchanID)]);
        
     
    %Calculate Theta ratio for plot/return    
    [thFFTspec,thFFTfreqs,t_FFT] = spectrogram(single(allLFP.data(:,goodTHidx)),window*Fs,noverlap*Fs,thFFTfreqs,Fs);
    t_FFT = t_FFT+allLFP.timestamps(1); %Offset for scoretime start
    t_spec = t_FFT;
    specdt = mode(diff(t_FFT));
    thFFTspec = (abs(thFFTspec));
    [zFFTspec,mu,sig] = zscore(log10(thFFTspec)');

    thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
    thpower = sum((thFFTspec(thfreqs,:)),1);
    allpower = sum((thFFTspec),1);

    thratio = thpower./allpower;    %Narrowband Theta
    thratio = smooth(thratio,smoothfact./specdt);
    
    %Remove ignoretimes (after smoothing), before normalizoing
    if ~isempty(ignoretime)
        ignoretimeIDX = InIntervals(t_FFT,ignoretime);
        thratio(ignoretimeIDX) = [];
        t_FFT(ignoretimeIDX) = [];
    end
    
subplot(5,1,3)
 %   plot(allLFP(:,1),allLFP(:,goodSWidx),'k')
    
        imagesc(t_spec,log2(thFFTfreqs),log10(thFFTspec))
        hold on
        plot(t_spec([1,end]),log2(f_theta([1,1])),'w')
        plot(t_spec([1,end]),log2(f_theta([2,2])),'w')
        axis xy
        plot(t_FFT,bz_NormToRange(thratio,log2(thFFTfreqs([1 end]))),'k','Linewidth',0.1)
        LogScale_ss('y',2)
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        ylim([log2(thFFTfreqs(1)) log2(thFFTfreqs(end))+0.2])
        ylabel({'LFP - FFT','f (Hz)'})
        title(['Theta Channel: ',num2str(THchanID)]);
        xlim(t_FFT([1,end]))
        set(gca,'XTick',[]);
        
        
        if saveFiles
saveas(chanfig,[figfolder,recordingname,'_SWTHChannels'],'jpeg')
        end
%saveas(chanfig,[figfolder,recordingname,'_SWTHChannels'],'fig')
end

