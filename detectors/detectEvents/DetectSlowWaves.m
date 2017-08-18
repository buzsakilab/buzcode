function [ SlowWaves,VerboseOut ] = DetectSlowWaves( basePath,varargin)
%[SlowWaves] = DetectSlowWaves(basePath) detects neocortical slow
%waves using a combination of a positive deflection in the LFP (delta wave)
%and a dip in gamma power.
%
%INPUTS
%   basePath  
%   (input options not yet implemented... all set to default)
%   'NREMInts'  -Interval of times for NREM 
%               -(Default: loaded from SleepState.states.mat)
%   'SWChan'   -Channel with the most robust (positively deflecting) Slow
%               Waves. (0-Indexing a la neuroscope). 
%               can try 'autoselect'
%               'useold' to use channel from existing SlowWave.states.mat
%               (Default: use SWChannel from sleep scoring)
%   'CTXSpkGroups' -Spike groups that are in the cortex...  default: all
%   'CTXChans' -LFP channels that are in the cortex...  default: all
%   'sensitivity' -sensititivity (0-1) for determining LFP thresholds
%                   gamma/delta thresholds set at the minimal peak 
%                   magnitude for which mean-normalized pop rate drops
%                   below the sensitivity value
%                   lower sensitivity will result in fewer False Positives,
%                   but more Missed slow waves. (default 0.6)
%   'filterparms' structure with fields
%           .deltafilter    [low high] bounds (default: [0.5 8]Hz)
%           .gammafilter    [low high] bounds (default: [100 400]Hz)
%           .gammasmoothwin  window for smoothing gamma power (default: 0.08 s)
%           .gammanormwin    window for normalizing gamma power (default: 20s)
%   'showFig'   -true/false show a quality control figure (default: true)
%   'saveMat'   -logical (default=true) to save in buzcode format
%   'forceReload' -logical (default: false) to redetect (add option to use
%                   old parameters like channels...)
%   'noPrompts'    -true/false disable any user prompts (default false)
%
%detectionparms : minOFF,mininterOFF, peakthresh (STD), peakdist (currently
%hardcoded)
%
%OUTPUTS
%   SlowWaves    a buzcode structure
%   VerboseOut   extra output stuff for detection quality checks/figures
%
%
%
%DLevenstein 2016/2017
%TO DO
%-incorporate multiple channels for detection of slow wave, which is robust
%on all (deep) lfp channels in the local cortical population
%update input parameter list
%% Defaults and Parms
ratevalidation = @(x) x>0 & x<1;
filterparmsvalidate = @(x) isstruct(x) & all(isfield(x,...
    {'deltafilter','gammafilter','gammasmoothwin','gammanormwin'}));

filterparms.deltafilter = [0.5 8];%heuristically defined.  room for improvement here.
filterparms.gammafilter = [100 400]; %high pass >80Hz (previously (>100Hz)
filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
filterparms.gammanormwin = 20; %window for gamma normalization (s)

p = inputParser;
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'showFig',true,@islogical);
addParameter(p,'SWChan',[]);
addParameter(p,'NREMInts',[]);
addParameter(p,'CTXChans','all');
addParameter(p,'sensitivity',0.6,ratevalidation);
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'filterparms',filterparms,filterparmsvalidate);
parse(p,varargin{:})

FORCEREDETECT = p.Results.forceReload;
SAVEMAT = p.Results.saveMat;
SHOWFIG = p.Results.showFig;
SWChann = p.Results.SWChan;
NREMInts = p.Results.NREMInts;
CTXChans = p.Results.CTXChans;
NOPROMPTS = p.Results.noPrompts;
ratethresh = p.Results.sensitivity;
filterparms = p.Results.filterparms;

%Defaults
if ~exist('basePath','var')
    basePath = pwd;
end

%Put this as optional input... with option to only use delta/gamma, for
%comparing quality


minwindur = 0.04;
joinwindur = 0.01;
%% File Management
baseName = bz_BasenameFromBasepath(basePath);
figfolder = fullfile(basePath,'DetectionFigures');
savefile = fullfile(basePath,[baseName,'.SlowWaves.events.mat']);
if exist(savefile,'file') && ~FORCEREDETECT
    display(['Slow Oscillation already Detected, loading ',baseName,'.SlowWaves.events.mat'])
    SlowWaves = bz_LoadEvents(basePath,'SlowWaves');
    return
end
%% Collect all the Necessary Pieces of Information
%Spikes in the CTX Spike Groups - assumes region ('CTX')
spikes = bz_GetSpikes('basepath',basePath,'region','CTX');
allspikes = sort(cat(1,spikes.times{:}));

%Sleep Scoring (for NREM)
[SleepState] = bz_LoadStates(basePath,'SleepState');
if isempty(SleepState) && isempty(NREMInts)
    button = questdlg(['SleepState.states.mat does not exist, '...
        'would you like to run SleepScoreMaster?'],...
        'DetectSlowWaves Needs NREM',...
        'Yes','No, use all timepoints','Cancel','Yes');
    switch button
        case 'Yes'
            SleepState = SleepScoreMaster(basePath);
            display('Please double check quality of sleep scoring in the StateScoreFigures folder')
            NREMInts = SleepState.ints.NREMstate;
        case 'No, use all timepoints'
            NREMInts = [0 Inf];
            SleepState.detectorparams.empty = [];
        case 'Cancel'
            return
    end
else
    NREMInts = SleepState.ints.NREMstate;
end
   
%LFP in the detection channel
if ~isfield(SleepState.detectorparams,'SWchannum') && isempty(SWChann)
    CHANSELECT = 'userinput';
    SWChann = input('Which channel shows the most robust (positive polarity) slow waves?');
elseif isempty(SWChann) %&& ismember(SleepState.detectorparams.SWchannum,CTXChans)
    CHANSELECT = 'SleepScoreSWchan';
    SWChann = SleepState.detectorparams.SWchannum;
elseif strcmp(SWChann,'manualselect')
    CHANSELECT = 'manual';
    %run SW channel selection routine: subfunction below (ManChanSelect)
elseif strcmp(SWChann,'autoselect')
    CHANSELECT = 'auto';
    %run SW channel selection routine: subfunction below (AutoChanSelect)
    SWChann = AutoChanSelect(CTXChans,basePath,NREMInts,spikes);
elseif strcmp(SWChann,'useold')
    SlowWaves = bz_LoadEvents(basePath,'SlowWaves');
    CHANSELECT = SlowWaves.detectorinfo.detectionparms.CHANSELECT;
    try
        SWChann = SlowWaves.detectorinfo.detectionchannel;
    catch %remove this tomorrow (8/15)
        SWChann = SlowWaves.detectorinfo.detectionparms.SWchannel;
    end
    clear SlowWaves
else
    CHANSELECT = 'userinput';
end
lfp = bz_GetLFP(SWChann,'basepath',basePath);

%% Filter the LFP: delta, gamma
display('Filtering LFP')
deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);
deltaLFP.normamp = NormToInt(deltaLFP.data,'modZ',NREMInts,deltaLFP.samplingRate);

gammaLFP = bz_Filter(lfp,'passband',filterparms.gammafilter,'filter','fir1','order',4);
gammaLFP.smoothamp = smooth(gammaLFP.amp,round(filterparms.gammasmoothwin.*gammaLFP.samplingRate),'moving' );
gammaLFP.normamp = NormToInt(gammaLFP.smoothamp,'modZ',NREMInts,gammaLFP.samplingRate,'moving',filterparms.gammanormwin);

%% Determine Thresholds and find windows around delta peaks/gamma dips
[thresholds,threshfigs] = DetermineThresholds(deltaLFP,gammaLFP,spikes,NREMInts,ratethresh,basePath);

%Find peaks in the delta LFP
[DELTApeaks,DELTAwins,DELTApeakheight] = FindPeakInWin(deltaLFP.normamp,deltaLFP.timestamps,...
    thresholds.DELTApeakthresh,thresholds.DELTAwinthresh,minwindur,joinwindur);
[DELTApeaks,keepPeaks] = RestrictInts(DELTApeaks,NREMInts);
DELTApeakheight = DELTApeakheight(keepPeaks);  DELTAwins = DELTAwins(keepPeaks,:);

%Find dips in the gamma power
[GAMMAdips,GAMMAwins,GAMMAdipdepth] = FindPeakInWin(-gammaLFP.normamp,gammaLFP.timestamps,...
    thresholds.GAMMAdipthresh,thresholds.GAMMAwinthresh,minwindur,joinwindur);
[GAMMAdips,keepPeaks] = RestrictInts(GAMMAdips,NREMInts);
GAMMAdipdepth = GAMMAdipdepth(keepPeaks);  GAMMAwins = GAMMAwins(keepPeaks,:);

%% Merge gamma/delta windows to get Slow Waves, UP/DOWN states
[ DOWNints,mergedidx ] = MergeSeparatedInts( [DELTAwins;GAMMAwins]);

%Keep only those windows in which DELTA/GAMMA were merged togehter (note
%this only works if joining happened previously... otherwise could keep
%windows where two delta and/or two gamma were joined; FindPeakInWin joins
numwins = cellfun(@length,mergedidx);
DOWNints = DOWNints(numwins>=2,:);
%Get the SW peak magnitude
mergedidx = mergedidx(numwins>=2); %keep only the indices from those that were merged into SWs
mergeddeltaidx = cellfun(@(X) X(X<=length(DELTAwins)),mergedidx,'UniformOutput',false); %keep delta
[SWpeakmag,peakidx] = cellfun(@(X) max(DELTApeakheight(X)),mergeddeltaidx,'UniformOutput',false); %Pick the larger of the peaks in each SW
SWpeaks = cellfun(@(X,Y) DELTApeaks(X(Y)),mergeddeltaidx,peakidx);
SWpeakmag = [SWpeakmag{:}]';

%% Remove UP/DOWN states that don't fit criteria

%Spiking in DOWN states
[status,interval,index] = InIntervals(allspikes,DOWNints); %spikes in DOWN
DOWNdur = diff(DOWNints,[],2);
DOWNrate = zeros(size(DOWNints,1),1);
for dd=1:size(DOWNints,1)
    DOWNrate(dd) = sum(interval==dd)./length(spikes.UID)./DOWNdur(dd);
end
aboveratethresh = DOWNrate>thresholds.ratethresh_Hz;
SWspikerejected = SWpeaks(aboveratethresh);
SWpeaks = SWpeaks(~aboveratethresh);
SWpeakmag = SWpeakmag(~aboveratethresh);
DOWNints = DOWNints(~aboveratethresh,:);
%Possible Issue here with onset/offset - if on/offset too late/early... 
%getting those spikes into the DOWN state

%Merge close DOWNs, take larger of the two peaks
[DOWNints,mergedidx] = MergeSeparatedInts(DOWNints,minwindur);
[SWpeakmag,newSWidx] = cellfun(@(X) max(SWpeakmag(X)),mergedidx,'UniformOutput',false);
newSWidx = cellfun(@(X,Y) X(Y),mergedidx,newSWidx);
SWpeaks = SWpeaks(newSWidx);
SWpeakmag = cat(1,SWpeakmag{:});

%Calculate UPs
UPints = [DOWNints(1:end-1,2) DOWNints(2:end,1)];
%Remove UPs that aren't fully in detectionints
UPints = RestrictInts(UPints,NREMInts);
UPdur = diff(UPints,[],2);


%%   Stuff for the detection figure
if SHOWFIG || nargout>1
%% Calculate Binned Population Rate
display('Binning Spikes')
numcells = length(spikes.UID);
dt = 0.005; %dt = 5ms
overlap = 8; %Overlap = 8 dt
winsize = dt*overlap; %meaning windows are 40ms big (previously 30)
[spikemat,t_spkmat,spindices] = SpktToSpkmat(spikes.times, [], dt,overlap);
synchmat = sum(spikemat>0,2);
ratemat = sum(spikemat,2);

%NREM spike rate histogram
[tspike_NREM,tidx_NREM] = RestrictInts(t_spkmat,NREMInts);
tspike_NREM = tspike_NREM(:,1);
synchmat_NREM = synchmat(tidx_NREM);
ratemat_NREM = ratemat(tidx_NREM);

%% CCG of SLow Waves and spikes
CCGvec = [SWpeaks,ones(size(SWpeaks));...
    allspikes,2.*ones(size(allspikes))];
[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');

%% Duration Histogram
UPDOWNdurhist.logbins = linspace(log10(0.02),log10(50),40);
[UPDOWNdurhist.DOWN,~] = hist(log10(DOWNdur),UPDOWNdurhist.logbins);
[UPDOWNdurhist.UP,~] = hist(log10(UPdur),UPDOWNdurhist.logbins);
UPDOWNdurhist.DOWN = UPDOWNdurhist.DOWN./sum(UPDOWNdurhist.DOWN);
UPDOWNdurhist.UP = UPDOWNdurhist.UP./sum(UPDOWNdurhist.UP);
%% Output for figues etc
VerboseOut.CCG.t_CCG =t_CCG;
VerboseOut.CCG.spikes = SWspike_CCG(:,1,2)./numcells;
VerboseOut.CCG.SWs = SWspike_CCG(:,1,1);
VerboseOut.UPDOWNdurhist = UPDOWNdurhist;
end

%% DETECTION FIGURE
if SHOWFIG
winsize = 8; %s
samplewin =randsample(SWpeaks,1)+winsize.*[-0.5 0.5];
sampleIDX = lfp.timestamps>=samplewin(1) & lfp.timestamps<=samplewin(2);

ratecolor = makeColorMap([1 1 1],[0.8 0 0],[0 0 0]);

figure('name',[baseName,' Slow Waves'])
    gammafig = copyobj(threshfigs.GAMMArate,gcf);
    subplot(4,2,4,gammafig)
    colorbar
    deltafig = copyobj(threshfigs.DELTArate,gcf);
    subplot(4,2,2,deltafig)
    colorbar
 
subplot(8,2,1)
    bar(t_CCG,SWspike_CCG(:,1,1),'facecolor','g','FaceAlpha',0.2)
    hold on
    plot([0 0],get(gca,'ylim'),'k-')
    xlim([-0.8 0.8])
    set(gca,'xticklabel',[]);
    ylabel('Delta Peaks')
    
subplot(8,2,3)
    bar(t_CCG,SWspike_CCG(:,1,2)./length(SWpeaks),'facecolor',[0.5 0.5 0.5])
    xlim([-0.8 0.8])
    set(gca,'xticklabel',[]);
    ylabel('Spike Rate');
    
subplot(4,2,3)
    plot(UPDOWNdurhist.logbins,UPDOWNdurhist.DOWN,'b','linewidth',2)
    hold on
    plot(UPDOWNdurhist.logbins,UPDOWNdurhist.UP,'r','linewidth',2)
    LogScale('x',10)
    axis tight
    legend('DOWN','UP')
    xlabel('Duration (s)') 
    
subplot(6,1,5)
    plot(lfp.timestamps(sampleIDX),lfp.data(sampleIDX),'k')
    axis tight
    hold on
    box off
    plot(lfp.timestamps(ismember(lfp.timestamps,SWpeaks)),lfp.data(ismember(lfp.timestamps,SWpeaks)),'g.')
    
        %DOWN patches
        y = get(gca,'Ylim');
        patch([DOWNints(:,1) DOWNints(:,2) DOWNints(:,2) DOWNints(:,1)]',...
            repmat([y(1) y(1) y(2) y(2)],length(DOWNints(:,1)),1)',...
            'g','FaceAlpha',0.2,'EdgeColor','none');
    
    xlim(samplewin)
    set(gca,'xticklabel',[]);  set(gca,'ytick',[])
    ylabel({'LFP',['Chan: ',num2str(SWChann)]})
    
subplot(6,1,4)
    bar(tspike_NREM,synchmat_NREM,'facecolor',[0.5 0.5 0.5])
    hold on
    %ylim([-8 length(unique(spikes.spindices(:,2)))+1])
    ylim([0 length(unique(spikes.spindices(:,2)))+1])
        %DOWN patches
        patch([DOWNints(:,1) DOWNints(:,2) DOWNints(:,2) DOWNints(:,1)]',...
            repmat([y(1) y(1) y(2) y(2)],length(DOWNints(:,1)),1)',...
            'g','FaceAlpha',0.2,'EdgeColor','none');
     
   % plot(deltaLFP.timestamps,deltaLFP.normamp-peakthresh,'b','Linewidth',1)    
   % plot(deltaLFP.timestamps(sampleIDX),gammaLFP.normamp(sampleIDX)+thresholds.GAMMAwinthresh,'g','Linewidth',1) 
   % plot(deltaLFP.timestamps(sampleIDX),deltaLFP.normamp(sampleIDX)-thresholds.DELTAwinthresh,'b','Linewidth',1) 
    plot(spikes.spindices(:,1),spikes.spindices(:,2),'k.','MarkerSize',4)
    plot(SWspikerejected,ones(size(SWspikerejected)),'kx','MarkerSize',10)
    box off    
    xlim(samplewin)
    ylabel('Cells');xlabel('t (s)')
    %legend('Spike Synchrony','Slow Waves','Gamma Power','Delta-Filtered LFP')
  
 subplot(6,1,6)
   % plot(deltaLFP.timestamps,deltaLFP.normamp-peakthresh,'b','Linewidth',1)    
   plot(deltaLFP.timestamps(sampleIDX),gammaLFP.normamp(sampleIDX),'g','Linewidth',1) 
   hold on
   plot(deltaLFP.timestamps(sampleIDX),deltaLFP.normamp(sampleIDX),'b','Linewidth',1) 
   axis tight
   box off
    xlim(samplewin)
    plot(GAMMAdips,-GAMMAdipdepth,'go')
plot(DELTApeaks,DELTApeakheight,'bo')
%plot(SWpeaks,SWpeakmag,'ko')
plot(GAMMAwins',-thresholds.GAMMAwinthresh.*ones(size(GAMMAwins')),'g')
plot(DELTAwins',thresholds.DELTAwinthresh.*ones(size(DELTAwins')),'b')

plot(samplewin,-thresholds.GAMMAwinthresh.*[1 1],'g:')
plot(samplewin,thresholds.DELTAwinthresh.*[1 1],'b:')
plot(samplewin,-thresholds.GAMMAdipthresh.*[1 1],'g--')
plot(samplewin,thresholds.DELTApeakthresh.*[1 1],'b--')

    
NiceSave('SlowOscillation',figfolder,baseName)
end

%% Ouput in .event.mat format
detectionparms.CHANSELECT = CHANSELECT;
detectionparms.CTXChans = CTXChans;
detectionparms.thresholds = thresholds;
detectionparms.filterparms = filterparms;

SlowWaves.timestamps = SWpeaks;
SlowWaves.SWpeakmag = SWpeakmag;
SlowWaves.ints.UP = UPints;
SlowWaves.ints.DOWN = DOWNints;
SlowWaves.detectorinfo.detectorname = 'DetectSlowWaves';
SlowWaves.detectorinfo.detectionparms = detectionparms;
SlowWaves.detectorinfo.detectiondate = today;
SlowWaves.detectorinfo.detectionintervals = NREMInts;
SlowWaves.detectorinfo.detectionchannel = SWChann;

if SAVEMAT
    save(savefile,'SlowWaves')
end

display('Slow Wave Detection: COMPLETE!')

if ~NOPROMPTS
    button = questdlg('Would you like to open EventExplorer to check Detection Quality?');
    switch button
        case 'Yes'
            EventExplorer(basePath,SlowWaves );
    end
end
%Prompt user here to check detection with EventExplorer
end








%% Channel Selection Functions
function usechan = AutoChanSelect(trychans,basePath,NREMInts,spikes)
    display('Detecting best channel for slow wave detection...')
    baseName = bz_BasenameFromBasepath(basePath);
    figfolder = fullfile(basePath,'DetectionFigures');
    %Exclude badchannels
    par = bz_getSessionInfo(basePath);
    if strcmp(trychans,'all') 
        trychans = [par.SpkGrps(:).Channels];
    end
    if isfield(par,'badchannels')
    	trychans = setdiff(trychans,par.badchannels);
    end
    
    %Calculate binned spike rate for correlation with gamma/anticorrelation
    %with LFP
    dt = 0.005; %dt = 5ms
    overlap = 8; %Overlap = 8 dt
    winsize = dt*overlap; %meaning windows are 40ms big (previously 30)
    [spikemat,t_spkmat,spindices] = SpktToSpkmat(spikes.times, [], dt,overlap);
    synchmat = sum(spikemat>0,2);
    ratemat = sum(spikemat,2);
    [t_spkmat,inNREMidx] = RestrictInts(t_spkmat,NREMInts); %Replace with InInterval
    synchmat = synchmat(inNREMidx);
    
  %%  
    for cc = 1:length(trychans)

        display(['Trying Channel ',num2str(cc),' of ',num2str(length(trychans))])
        %Load the LFPs
        chanlfp = bz_GetLFP(trychans(cc),'basepath',basePath);
        %Filter in gamma
        gammafilter = [100 512];
        trygammaLFP = bz_Filter(chanlfp,'passband',gammafilter,'order',4);

        %Restrict to NREM only - could also use intervals above to do this....
        [chanlfp.timestamps,inNREMidx] = RestrictInts(chanlfp.timestamps,NREMInts);
        chanlfp.data = chanlfp.data(inNREMidx,:);
        trygammaLFP.amp = trygammaLFP.amp(inNREMidx,:);
        %Best channel is the one in which gamma is most anticorrelated with the
        %lowpass LFP - i.e. DOWN states (positive LFP) have low gamma power
        lowpassLFP = bz_Filter(chanlfp,'passband',[0 6],'order',2);
        
        gammaLFPcorr(cc) = corr(lowpassLFP.data,trygammaLFP.amp,'type','spearman');
        
        %Find LFP at rate time points and calculate correlation
        lowpasslfp.ratetimes = interp1(lowpassLFP.timestamps,lowpassLFP.data,t_spkmat,'nearest');
        trygammaLFP.ratetimes = interp1(lowpassLFP.timestamps,trygammaLFP.amp,t_spkmat,'nearest');
        lowpassspikecorr(cc) = corr(lowpasslfp.ratetimes,synchmat,'type','spearman');
        gammaspikecorr(cc) = corr(trygammaLFP.ratetimes,synchmat,'type','spearman');
        
       
        %Save a small window for plotting
        if ~exist('samplewin','var')
            winsize = 4; %s
            samplewin =randsample(chanlfp.timestamps,1)+winsize.*[-0.5 0.5];
            sampleIDX = chanlfp.timestamps>=samplewin(1) & chanlfp.timestamps<=samplewin(2);
            alllfp.timestamps = chanlfp.timestamps(sampleIDX);
        end
        alllfp.data(:,cc) = chanlfp.data(sampleIDX);
        alllfp.channels(cc) = chanlfp.channels;
        
    end

    
    %%
    %[~,usechanIDX] = min(gammaLFPcorr);
    [~,usechanIDX] = min(lowpassspikecorr.*gammaspikecorr);
    usechan = trychans(usechanIDX);
    
    display(['Selected Channel: ',num2str(usechan)])
    
    %[~,sortcorr] = sort(gammaLFPcorr);
    [~,sortcorr] = sort(lowpassspikecorr.*gammaspikecorr);
    %% Figure

    
    figure('name',[baseName,' Slow Wave Channel Selection'])
    subplot(2,2,1)
        hist(gammaLFPcorr)
        xlabel('Gamma-LFP Correlation')
%     subplot(2,2,3)
%         plot(trygammaLFP.amp(:,usechanIDX),chanlfp.data(:,usechanIDX),'k.')
%         xlabel('Gamma Power');ylabel('Raw LFP')
    subplot(6,2,4:2:12)
        bz_MultiLFPPlot(alllfp,'channels',trychans(sortcorr),'timewin',samplewin)
        ylabel('Slow Wave Channel (Worse<--------->Better)')
        
    subplot(6,2,2)
        plot(alllfp.timestamps,alllfp.data(:,usechanIDX),'k')
        axis tight
        box off
        
    subplot(2,2,3)
    plot(lowpassspikecorr,gammaspikecorr,'k.')
    hold on
    	plot(lowpassspikecorr(usechanIDX),gammaspikecorr(usechanIDX),'ro')
        xlabel('Spike-Delta Correlation');ylabel('Spike-Gamma Power Correlation')
        
    NiceSave('SlowWaveChannelSelect',figfolder,baseName)
end



%% Threshold Determination Function
function [thresholds,threshfigs] = DetermineThresholds(deltaLFP,gammaLFP,spikes,NREMInts,ratethresh,basePath)
    display('Determining delta/gamma thresholds for detection...')
    baseName = bz_BasenameFromBasepath(basePath);
    minwindur = 0.04; %should pass through... but not super important
    
    %Find peaks in delta, gamma power
    [peakheights,DELTApeaks] = findpeaks(deltaLFP.normamp,deltaLFP.timestamps,'MinPeakHeight',0.25,'MinPeakDistance',minwindur);
    [DELTApeaks,keepPeaks] = RestrictInts(DELTApeaks,NREMInts);
    DELTAPeakheight = peakheights(keepPeaks);
    [~,DELTApeakIDX] = ismember(DELTApeaks,deltaLFP.timestamps);

    [peakheights,GAMMAdips] = findpeaks(-gammaLFP.normamp,gammaLFP.timestamps,'MinPeakHeight',0.5,'MinPeakDistance',minwindur);
    [GAMMAdips,keepPeaks] = RestrictInts(GAMMAdips,NREMInts);
    GAMMAdipdepth = peakheights(keepPeaks);
    [~,GAMMAdipIDX] = ismember(GAMMAdips,gammaLFP.timestamps);
    
    %% Get the spike PETH around delta/gamma peaks
    display('Calculating PETH by Peak For Threshold Calibration')
    win = 1;
    numchecks = 20000;

    numtimebins = 200;
    nummagbins = 30;
    timebins = linspace(-win,win,numtimebins);
    
    numcells = length(spikes.times);
    allspikes = sort(cat(1,spikes.times{:}));
    
    %DELTA
    display('DELTA...')
    reltime = [];
    spkpeakheight = [];
    if length(DELTApeaks)>numchecks
        sampleDELTA = randsample(length(DELTApeaks),numchecks);
    else
        sampleDELTA = 1:length(DELTApeaks);
    end
    nearpeakdelta = zeros(2.*win.*deltaLFP.samplingRate+1,nummagbins);
    DELTAmagbins = linspace(min(DELTAPeakheight),min([max(DELTAPeakheight)-1,5]), nummagbins);
    for pp = 1:length(sampleDELTA)
        ss = sampleDELTA(pp);
        %Spikes around the delta
        nearpeakspikes = allspikes >= DELTApeaks(ss)-win & allspikes <= DELTApeaks(ss)+win;
        reltime = [reltime; allspikes(nearpeakspikes)-DELTApeaks(ss)];
        spkpeakheight = [spkpeakheight; DELTAPeakheight(ss).*ones(sum(nearpeakspikes),1)];
        %Delta around the delta
        [~,groupidx] = min(abs(DELTAPeakheight(ss)-DELTAmagbins));
        nearpeakdelta(:,groupidx) =nansum([nearpeakdelta(:,groupidx), ...
            deltaLFP.normamp(DELTApeakIDX(ss)+[-1.*deltaLFP.samplingRate:deltaLFP.samplingRate])],2);
    end
    peakmagdist = hist(DELTAPeakheight(sampleDELTA),DELTAmagbins);
    spikehitmat = hist3([reltime,spkpeakheight],{timebins,DELTAmagbins});
    ratemat_byDELTAmag = bsxfun(@(A,B) A./B./mean(diff(timebins))./numcells,spikehitmat,peakmagdist);
    deltapower_byDELTAmag = bsxfun(@(A,B) A./B,nearpeakdelta,peakmagdist);
    
    %GAMMA
    display('GAMMA...')
    reltime = [];
    spkpeakheight = [];
    if length(GAMMAdips)>numchecks
        sampleGAMMA = randsample(length(GAMMAdips),numchecks); %Should really try to get even# in each bin...
    else
        sampleGAMMA = 1:length(GAMMAdips);
    end
    GAMMAmagbins = linspace(min(GAMMAdipdepth),min([max(GAMMAdipdepth)-0.5,2.5]),nummagbins);
    neardipgamma = zeros(2.*win.*gammaLFP.samplingRate+1,nummagbins);
    for pp = 1:length(sampleGAMMA)
        ss = sampleGAMMA(pp);
        %Spikes around the gamma
        nearpeakspikes = allspikes >= GAMMAdips(ss)-win & allspikes <= GAMMAdips(ss)+win;
        reltime = [reltime; allspikes(nearpeakspikes)-GAMMAdips(ss)];
        spkpeakheight = [spkpeakheight; GAMMAdipdepth(ss).*ones(sum(nearpeakspikes),1)];
         %Gamma around the gamma
        [~,groupidx] = min(abs(GAMMAdipdepth(ss)-GAMMAmagbins)); %Which gamma power bin
        neardipgamma(:,groupidx) =nansum([neardipgamma(:,groupidx), ... %nans for non-NREM from window smoothing
            gammaLFP.normamp(GAMMAdipIDX(ss)+[-1.*gammaLFP.samplingRate:gammaLFP.samplingRate])],2);
    end
    peakmagdist = hist(GAMMAdipdepth(sampleGAMMA),GAMMAmagbins);
    spikehitmat = hist3([reltime,spkpeakheight],{timebins,GAMMAmagbins});
    ratemat_byGAMMAmag = bsxfun(@(A,B) A./B./mean(diff(timebins))./numcells,spikehitmat,peakmagdist);
    gammapower_byGAMMAmag = bsxfun(@(A,B) A./B,neardipgamma,peakmagdist);
    %% Find Thresholds as values where spike rate falls below threshold
    %Set the rate threshold to halfway between that at the delta peak and
    %the average around DELTA
    ratemat_byDELTAmag = imgaussfilt(ratemat_byDELTAmag,2);
    ratemat_byGAMMAmag = imgaussfilt(ratemat_byGAMMAmag,2);
    meanratearoundDELTA = mean(ratemat_byDELTAmag(:));
    minrateatDELTApeak = min(ratemat_byDELTAmag(round(end/2),:));
    DELTAraterange = meanratearoundDELTA-minrateatDELTApeak;
    
    ratethresh_Hz = minrateatDELTApeak+(ratethresh.*DELTAraterange);

    DELTAbox=bwmorph(ratemat_byDELTAmag<ratethresh_Hz,'close');
    DELTAbox=bwmorph(DELTAbox,'open');
    if sum(DELTAbox(:))== 0
        display('No DOWN around gamma dip.... perhaps adjust rate threshold or pick another channel?')
        DELTApeakthresh = 2.2;
        DELTAwinthresh = 1;
        DELTAbox = [1 1];
    else
    DELTAbox=bwboundaries(DELTAbox); DELTAbox = DELTAbox{1};
    DELTAbox(DELTAbox(:,2)==max(DELTAbox(:,2)),:)=[];
    DELTApeakthresh = DELTAmagbins(min(DELTAbox(:,2)));
    scalefactor = length(nearpeakdelta)./numtimebins; %convert number of bins for LFP and spikes
    DELTAwinthresh = deltapower_byDELTAmag(sub2ind(size(deltapower_byDELTAmag),...
        round(DELTAbox(:,1).*scalefactor),DELTAbox(:,2)));
    DELTAwinthresh = mean(DELTAwinthresh);
    end
    
    GAMMAbox=bwmorph(ratemat_byGAMMAmag<ratethresh_Hz,'close');
    GAMMAbox=bwmorph(GAMMAbox,'open');
    if sum(GAMMAbox(:))== 0   %will have issue here with no dip recordings... bad channel.
        display('No DOWN around gamma dip.... perhaps adjust rate threshold or pick another channel?')
        GAMMAdipthresh = 1.2;
        GAMMAwinthresh = 1;
        GAMMAbox = [1 1];
    else
    GAMMAbox=bwboundaries(GAMMAbox); GAMMAbox = GAMMAbox{1}; %might have issue here with no dip recordings... "multiple objects"
    GAMMAbox(GAMMAbox(:,2)==max(GAMMAbox(:,2)),:)=[];
    GAMMAdipthresh = GAMMAmagbins(min(GAMMAbox(:,2)));
    scalefactor = length(neardipgamma)./numtimebins; %convert number of bins for LFP and spikes
    GAMMAwinthresh = gammapower_byGAMMAmag(sub2ind(size(gammapower_byGAMMAmag),...
        round(GAMMAbox(:,1).*scalefactor),GAMMAbox(:,2)));
    GAMMAwinthresh = -mean(GAMMAwinthresh);
    end
    
    thresholds.DELTApeakthresh = DELTApeakthresh;
    thresholds.DELTAwinthresh = DELTAwinthresh;
    thresholds.GAMMAdipthresh = GAMMAdipthresh;
    thresholds.GAMMAwinthresh = GAMMAwinthresh;
    thresholds.ratethresh = ratethresh;
    thresholds.ratethresh_Hz = ratethresh_Hz;

    %%
    %rate normalization for plots
    normrate_DELTA = (ratemat_byDELTAmag-minrateatDELTApeak)./DELTAraterange;
    normrate_GAMMA = (ratemat_byGAMMAmag-minrateatDELTApeak)./DELTAraterange;
    
    figure('name',[baseName,' Threshold Detection'])
     threshfigs.DELTArate = subplot(2,2,3);
        imagesc(timebins,DELTAmagbins,normrate_DELTA')
        hold on
        plot(timebins(DELTAbox(:,1)),DELTAmagbins(DELTAbox(:,2)),'r.')
        plot(get(gca,'xlim'),DELTApeakthresh.*[1 1],'k--')
        plot(get(gca,'xlim'),DELTAwinthresh.*[1 1],'k:')
        caxis([0 1.5])
        colorbar
        xlabel('t (relative to SW peak)');ylabel({'SW Peak Amplitude', '(modZ)'})
        axis xy
        xlim([-0.7 0.7])

    threshfigs.GAMMArate = subplot(2,2,4);
        imagesc(timebins,GAMMAmagbins,normrate_GAMMA')
        hold on
        plot(timebins(GAMMAbox(:,1)),GAMMAmagbins(GAMMAbox(:,2)),'r.')
        plot(get(gca,'xlim'),GAMMAdipthresh.*[1 1],'k--')
        plot(get(gca,'xlim'),GAMMAwinthresh.*[1 1],'k:')
        caxis([0 1.5])
        colorbar
        xlabel('t (relative to GA Dip)');ylabel({'GA Dip Amplitude', '(modZ)'})
        xlim([-0.7 0.7])

    threshfigs.DELTApower = subplot(2,2,1);
        imagesc(timebins,DELTAmagbins,deltapower_byDELTAmag')
        hold on
        plot(timebins(DELTAbox(:,1)),DELTAmagbins(DELTAbox(:,2)),'r.')
        plot(get(gca,'xlim'),DELTApeakthresh.*[1 1],'k--')
        plot(get(gca,'xlim'),DELTAwinthresh.*[1 1],'k:')
        colorbar
        xlabel('t (relative to SW peak)');ylabel({'SW Peak Amplitude', '(modZ)'})
        axis xy
        xlim([-0.7 0.7])

    threshfigs.GAMMApower = subplot(2,2,2);
        imagesc(timebins,GAMMAmagbins,gammapower_byGAMMAmag')
        hold on
        plot(timebins(GAMMAbox(:,1)),GAMMAmagbins(GAMMAbox(:,2)),'r.')
        plot(get(gca,'xlim'),GAMMAdipthresh.*[1 1],'k--')
        plot(get(gca,'xlim'),GAMMAwinthresh.*[1 1],'k:')
        xlabel('t (relative to GA Dip)');ylabel({'GA Dip Amplitude', '(modZ)'})
        colorbar
        xlim([-0.7 0.7])
        
end


