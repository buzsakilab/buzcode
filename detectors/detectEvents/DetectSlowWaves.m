function [ SlowWaves,VerboseOut ] = DetectSlowWaves( basePath,varargin)
%[SlowWaves] = DetectSlowWaves(basePath) detects neocortical slow
%waves using a combination of a positive deflection in the LFP (delta wave)
%and a dip in gamma power.
%
%INPUTS
%   basePath    -(default: pwd) basePath for the recording file, in buzcode format:
%                   whateverPath/baseName/
%               folder must include files:
%                   baseName.lfp
%                   baseName.sessionInfo.mat -or- baseName.xml
%               recommended (will prompt if not included unless 'noPrompts')
%                   baseName.spikes.cellinfo.mat (or clu/res/fet)
%                   baseName.SleepState.states.mat
%               
%   (options)
%   'lfp'               -A buzcode-style lfp structure... if you would
%                        rather just input the lfp instead of loading from
%                        basepath
%                           Default: load from basePath with bz_GetLFP
%   'spikes'            -A buzcode-style spike structure 
%                           Default: load from basePath with bz_GetSpikes
%   'NREMInts'          -Interval of times for NREM (seconds) 
%                        (Default: loaded from SleepState.states.mat, 
%                                   run SleepScoreMaster if not exist)
%                        use [0 Inf] to detect over all time points
%   'DetectionChannel'  -Channel with the most robust Slow Waves. (0-Indexing a la neuroscope). 
%                        (Default: 'autoselect')
%                        'useold' to use channel from existing SlowWaves.events.mat 
%                        If providing lfp via the 'lfp' input, make sure to
%                        give a detection channel here.
%   'noSpikes'          -true/false - set to true to not use spike information
%                        (default: false)
%   'MUAspikes'       -true/false - use MUA peaks (500-5000Hz) extracted 
%                        from the .dat file instead of spikes
%   'CTXChans'          -LFP channels that are in the cortex...  
%                        default: region 'CTX' from baseName.sessionInfo.mat or xml
%   'sensitivity'       -sensititivity (0-1) for determining LFP thresholds
%                        sensitivity for setting gamma/delta thresholds.
%                        lower sensitivity will result in fewer False Positives,
%                        but more Missed slow waves. (default 0.6)
%   'filterparms'       -filtering parameters structure with fields:
%           .deltafilter    [low high] bounds (default: [0.5 8]Hz)
%           .gammafilter    [low high] bounds (default: [100 400]Hz)
%           .gammasmoothwin  window for smoothing gamma power (default: 0.08 s)
%           .gammanormwin    window for normalizing gamma power (default: 20s)
%   'showFig'           -true/false show a quality control figure (default: true)
%   'saveMat'           -logical (default=true) to save in buzcode format
%   'forceReload'       -logical (default: false) to re-detect
%   'noPrompts'         -true/false disable any user prompts (default: false)
%
%
%OUTPUTS
%   SlowWaves    a buzcode structure
%   VerboseOut   extra output stuff for detection quality checks/figures
%
%
%
%DLevenstein 2016/2017
%If used, please cite: Levenstein et al 2018, currently on bioRxiv
%TO DO
%-incorporate multiple channels for detection of slow wave, which is robust
%on all (deep) lfp channels in the local cortical population
%-update input parameter list
%NOTE: requires 2017a or higher. sad. (functions: movingmad and movingmedian)
%% Defaults and Parms
ratevalidation = @(x) x>0 & x<1;
filterparmsvalidate = @(x) isstruct(x) & all(isfield(x,...
    {'deltafilter','gammafilter','gammasmoothwin','gammanormwin'}));

filterparms.deltafilter = [0.5 8];%heuristically defined.  room for improvement here.
filterparms.gammafilter = [100 400];
filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
filterparms.gammanormwin = 20; %window for gamma normalization (s)

p = inputParser;
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'showFig',true,@islogical);
addParameter(p,'noSpikes',false,@islogical);
addParameter(p,'MUAspikes',false,@islogical);
addParameter(p,'DetectionChannel','autoselect');
addParameter(p,'NREMInts',[]);
addParameter(p,'lfp',[]);
addParameter(p,'spikes',[]);
addParameter(p,'CTXChans','all');
addParameter(p,'sensitivity',0.6,ratevalidation);
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'filterparms',filterparms,filterparmsvalidate);
addParameter(p,'minwindur',0.04);
addParameter(p,'joinwindur',0.01);


parse(p,varargin{:})

FORCEREDETECT = p.Results.forceReload;
SAVEMAT = p.Results.saveMat;
SHOWFIG = p.Results.showFig;
SWChan = p.Results.DetectionChannel;
NREMInts = p.Results.NREMInts;
CTXChans = p.Results.CTXChans;
noPrompts = p.Results.noPrompts;
NOSPIKES = p.Results.noSpikes;
ratethresh = p.Results.sensitivity;
filterparms = p.Results.filterparms;
lfp = p.Results.lfp;
spikes = p.Results.spikes;
MUAspikes = p.Results.MUAspikes;
minwindur = p.Results.minwindur;
joinwindur = p.Results.joinwindur;

%Defaults
if ~exist('basePath','var')
    basePath = pwd;
end

%% File Management
baseName = bz_BasenameFromBasepath(basePath);
figfolder = fullfile(basePath,'DetectionFigures');
savefile = fullfile(basePath,[baseName,'.SlowWaves.events.mat']);
if exist(savefile,'file') && ~FORCEREDETECT
    display(['Slow Oscillation already Detected, loading ',baseName,'.SlowWaves.events.mat'])
    SlowWaves = bz_LoadEvents(basePath,'SlowWaves');
    return
end
%% Collect all the Necessary Pieces of Information: Spikes, States, LFP

%Spikes in the CTX Spike Groups - assumes region ('CTX')
if NOSPIKES
    spikes = 'NOSPIKES'; allspikes = 'NOSPIKES';
    numcells = nan;
elseif ~isempty(spikes)
    allspikes = sort(cat(1,spikes.times{:}));
    numcells = length(spikes.UID);
elseif MUAspikes
    [ MUA ] = MUAfromDat( basePath,'usepeaks',true,'saveMat',true );
    spikes = MUA.peaks;
    allspikes = sort(cat(1,spikes.times{:}));
    numcells = length(spikes.times);
else
    try spikes = bz_GetSpikes('basepath',basePath,'region','CTX','noPrompts',noPrompts);
    catch
        spikes = bz_GetSpikes('basepath',basePath,'noPrompts',noPrompts);
    end
    if isempty(spikes)
        button = questdlg({['No spikes found (baseName.spikes.cellinfo.mat or clu/res/fet), '...
            'would you like to run in ''noSpikes'' mode?'],...
            '(Not recommended, lower sensitivity neded for good detection quality with no spikes)'},...
            'DetectSlowWaves Needs Spikes',...
            'Yes','Yes, with lower sensitivity (0.4)','No','Yes');
        if strcmp(button,'Yes')
            spikes = 'NOSPIKES'; allspikes = 'NOSPIKES';
            NOSPIKES = true;
            numcells = nan;
        elseif strcmp(button,'Yes, with lower sensitivity (0.4)')
            spikes = 'NOSPIKES'; allspikes = 'NOSPIKES';
            NOSPIKES = true;
            ratethresh = 0.4;
            numcells = nan;
        end
    else
        allspikes = sort(cat(1,spikes.times{:}));
        numcells = length(spikes.UID);
    end

end

%Sleep Scoring (for NREM). Load, Prompt if doesn't exist
if isempty(NREMInts)
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
end
   
%Which channel should be used for detection?
switch SWChan
    case 'autoselect'	%Automated selection
        CHANSELECT = 'auto';
        [SWChan,CTXChans] = AutoChanSelect(CTXChans,basePath,NREMInts,spikes,filterparms,noPrompts);
    case 'useold'       %Use from existing SlowWaves.events.mat
        SlowWaves = bz_LoadEvents(basePath,'SlowWaves');
        CHANSELECT = SlowWaves.detectorinfo.detectionparms.CHANSELECT;
        SWChan = SlowWaves.detectorinfo.detectionchannel;
        clear SlowWaves
    otherwise 
        if isnumeric(SWChan)    %If user has entered channel number
            CHANSELECT = 'userinput';
        else
            SWChan = inputdlg(['Which LFP channel would you like to ',...
                'use for detection? (0-index a la neuroscope)'],...
                'DetectionChannel',1);
            CHANSELECT = 'userinput';
            assert(isnumeric(SWChan),'Your SW channel should be a number, eh?');
        end     %Add: option to use SleepState.detectorparams.SWchannum?
end
        
%Load the LFP
if isempty(lfp)
    lfp = bz_GetLFP(SWChan,'basepath',basePath,'noPrompts',noPrompts);
end

%% Filter the LFP: delta, high gamma 
display('Filtering LFP')
deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);
deltaLFP.normamp = NormToInt(deltaLFP.data,'modZ',NREMInts,deltaLFP.samplingRate);

%switch useMUA
%    case false
gammaLFP = bz_Filter(lfp,'passband',filterparms.gammafilter,'filter','fir1','order',4);
%    case true
%         [ MUA ] = MUAfromDat( basePath,'channels',SWChan);
%         gammaLFP.amp = MUA.data;
%         gammaLFP.samplingRate = MUA.samplingRate;
%         gammaLFP.timestamps = MUA.timestamps;
% end
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
[ DOWNints,mergedidx ] = MergeSeparatedInts([DELTAwins;GAMMAwins]);

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
DOWNdur = diff(DOWNints,[],2);

if ~NOSPIKES
%Spiking in DOWN states
    [status,interval,index] = InIntervals(allspikes,DOWNints); %spikes in DOWN
    DOWNrate = zeros(size(DOWNints,1),1);
    for dd=1:size(DOWNints,1)
        DOWNrate(dd) = sum(interval==dd)./numcells./DOWNdur(dd);
    end
    aboveratethresh = DOWNrate>thresholds.ratethresh_Hz;
    SWspikerejected = SWpeaks(aboveratethresh);
    SWpeaks = SWpeaks(~aboveratethresh);
    SWpeakmag = SWpeakmag(~aboveratethresh);
    DOWNints = DOWNints(~aboveratethresh,:);
    %Possible Issue here with onset/offset - if on/offset too late/early... 
    %getting those spikes into the DOWN state
end

%Merge close DOWNs, take larger of the two peaks
[DOWNints,mergedidx] = MergeSeparatedInts(DOWNints,joinwindur);
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
%% Duration Histogram
UPDOWNdurhist.logbins = linspace(log10(0.02),log10(50),40);
[UPDOWNdurhist.DOWN,~] = hist(log10(DOWNdur),UPDOWNdurhist.logbins);
[UPDOWNdurhist.UP,~] = hist(log10(UPdur),UPDOWNdurhist.logbins);
UPDOWNdurhist.DOWN = UPDOWNdurhist.DOWN./sum(UPDOWNdurhist.DOWN);
UPDOWNdurhist.UP = UPDOWNdurhist.UP./sum(UPDOWNdurhist.UP);
    
%% Calculate Binned Population Rate
switch NOSPIKES
    case false
display('Binning Spikes')
overlap = 6; %Overlap = 8 dt
binsize = 0.03; %meaning windows are 40ms big (previously 30)
spikemat = bz_SpktToSpkmat(spikes, 'binsize', binsize,'overlap',overlap);
synchmat = sum(spikemat.data>0,2);
ratemat = sum(spikemat.data,2);

%NREM spike rate histogram
tidx_NREM = InIntervals(spikemat.timestamps,NREMInts);
tspike_NREM = spikemat.timestamps(tidx_NREM);
synchmat_NREM = synchmat(tidx_NREM);
ratemat_NREM = ratemat(tidx_NREM);

%% CCG of SLow Waves and spikes
CCGvec = [SWpeaks,ones(size(SWpeaks));...
    allspikes,2.*ones(size(allspikes))];
[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');

    case true
tspike_NREM = nan;
synchmat_NREM = nan;
SWspikerejected = nan;
clear spikes
spikes.spindices = [nan nan];
        
CCGvec = [SWpeaks,ones(size(SWpeaks));...
    SWpeaks,2.*ones(size(SWpeaks));];
[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02,'norm','rate');


end

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
    ylabel({'LFP',['Chan: ',num2str(SWChan)]})
    
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

try    
    NiceSave('SlowOscillation',figfolder,baseName)
catch
    display('ERROR SAVING FIGURE')
end
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
SlowWaves.detectorinfo.detectiondate = datetime('today');
SlowWaves.detectorinfo.detectionintervals = NREMInts;
SlowWaves.detectorinfo.detectionchannel = SWChan;

try
    bz_tagChannel(basePath,SWChan,'SWChan','noPrompts',noPrompts);
catch
    display('Unable to save channel tag in sessionInfo')
end

if SAVEMAT
    save(savefile,'SlowWaves')
end

display('Slow Wave Detection: COMPLETE!')

if ~noPrompts
    button = questdlg('Would you like to open EventExplorer to check Detection Quality?');
    switch button
        case 'Yes'
            EventExplorer(basePath,SlowWaves);
    end
end

end








%% Channel Selection Functions
function [usechan,trychans] = AutoChanSelect(trychans,basePath,NREMInts,spikes,filterparms,noPrompts)
    display('Detecting best channel for slow wave detection...')
    baseName = bz_BasenameFromBasepath(basePath);
    figfolder = fullfile(basePath,'DetectionFigures');
    
    par = bz_getSessionInfo(basePath,'noPrompts',noPrompts); %Load the metadata
    if strcmp(trychans,'all') %ie the user hasn't manually entered channels
        try %In case there is no region field in the SessionInfo
            regions = unique(par.region);
        catch
            regions = {''};
        end
        if length(regions)==1 && strcmp(regions,'')
            %If there are no regions, use all (tell user)
            display('No region in sessionInfo, using all channels')
            trychans = par.channels;
        elseif length(regions)>0 && any(strcmp(regions,'CTX')) 
            %If there's a region 'CTX', use those channels
            trychans = par.channels(strcmp(par.region,'CTX'));
        elseif length(regions)>0 && ~any(strcmp(regions,'CTX')) 
            %If there are regions, but not CTX, allow user to selet regions
            [s,v] = listdlg('PromptString','No ''CTX'', which region to use? :',...
                'SelectionMode','multiple',...
                'ListString',regions)
            trychans = par.channels(ismember(par.region,regions(s)));
        end
    end
    if isfield(par,'badchannels')     %Exclude badchannels
    	trychans = setdiff(trychans,par.badchannels);
    end
    
    %Calculate binned spike rate for correlation with gamma/anticorrelation
    %with LFP
    if strcmp(spikes,'NOSPIKES'); NOSPIKES=true; else NOSPIKES=false; end
    if ~NOSPIKES
        dt = 0.005; %dt = 5ms
        overlap = 8; %Overlap = 8 dt
        winsize = dt*overlap; %meaning windows are 40ms big (previously 30)
        %[spikemat,t_spkmat,spindices] = bz_SpktToSpkmat(spikes.times, [], dt,overlap);
        spikemat = bz_SpktToSpkmat(spikes, 'binsize', winsize,'overlap',overlap);
        synchmat = sum(spikemat.data>0,2);
        ratemat = sum(spikemat.data,2);
        [t_spkmat,inNREMidx] = RestrictInts(spikemat.timestamps,NREMInts); %Replace with InInterval
        synchmat = synchmat(inNREMidx);
    end
    
  %%  
    for cc = 1:length(trychans)
        display(['Trying Channel ',num2str(cc),' of ',num2str(length(trychans))])
        chanlfp = bz_GetLFP(trychans(cc),'basepath',basePath,'noPrompts',true);
        %Filter in gamma
        gammafilter = filterparms.gammafilter; %Note: this doesn't work as well with new filtered LFP.... need better MUA
        trygammaLFP = bz_Filter(chanlfp,'passband',gammafilter,'order',4);

        %Restrict to NREM only - could also use intervals above to do this....
        [chanlfp.timestamps,inNREMidx] = RestrictInts(chanlfp.timestamps,NREMInts);
        chanlfp.data = chanlfp.data(inNREMidx,:);
        trygammaLFP.amp = trygammaLFP.amp(inNREMidx,:);
        %Best channel is the one in which gamma is most anticorrelated with the
        %lowpass LFP - i.e. DOWN states (positive LFP) have low gamma power
        lowpassLFP = bz_Filter(chanlfp,'passband',filterparms.deltafilter,'order',2);
        
        gammaLFPcorr(cc) = corr(lowpassLFP.data,trygammaLFP.amp,'type','spearman');
        
        %Find LFP at rate time points and calculate correlation
        if ~NOSPIKES
            lowpasslfp.ratetimes = interp1(lowpassLFP.timestamps,lowpassLFP.data,t_spkmat,'nearest');
            trygammaLFP.ratetimes = interp1(lowpassLFP.timestamps,trygammaLFP.amp,t_spkmat,'nearest');
            lowpassspikecorr(cc) = corr(lowpasslfp.ratetimes,synchmat,'type','spearman');
            gammaspikecorr(cc) = corr(trygammaLFP.ratetimes,synchmat,'type','spearman');
        end
        
       
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
    if NOSPIKES
        [~,usechanIDX] = min(gammaLFPcorr);
        [~,sortcorr] = sort(gammaLFPcorr);
        lowpassspikecorr = nan;gammaspikecorr = nan;
    else
        [~,usechanIDX] = min(lowpassspikecorr.*gammaspikecorr);
        [~,sortcorr] = sort(lowpassspikecorr.*gammaspikecorr);
    end
    usechan = trychans(usechanIDX);
    
    display(['Selected Channel: ',num2str(usechan)])
    
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
    
    if ~NOSPIKES  %the following plot relies on spikes  
    subplot(2,2,3)
        plot(lowpassspikecorr,gammaspikecorr,'k.')
        hold on
    	plot(lowpassspikecorr(usechanIDX),gammaspikecorr(usechanIDX),'ro')
        xlabel('Spike-Delta Correlation');ylabel('Spike-Gamma Power Correlation')
    end
        
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

    [peakheights,GAMMAdips] = findpeaks(-gammaLFP.normamp,gammaLFP.timestamps,'MinPeakHeight',0.2,'MinPeakDistance',minwindur);
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
    
    if strcmp(spikes,'NOSPIKES'); 
        NOSPIKES=true; allspikes = [];
    else
        numcells = length(spikes.times);
        allspikes = sort(cat(1,spikes.times{:}));
        NOSPIKES=false;
    end
    
    %DELTA
    display('DELTA...')
    %Pick the delta peaks to average spiking around
    if length(DELTApeaks)>numchecks
        sampleDELTA = randsample(length(DELTApeaks),numchecks);
    else
        sampleDELTA = 1:length(DELTApeaks);
    end
    %Remove sample deltas too close to begining/end of recording
        sampleDELTA(DELTApeaks(sampleDELTA)<=win) = [];
        sampleDELTA(DELTApeaks(sampleDELTA)>=(deltaLFP.timestamps(end)-win)) = [];
    %Loop Through delta peaks and calculate spike rate (or high gamma)
    %around them, as a function of their peak height
    relspktime = [];    %Spike time relative to delta peaks
    spkpeakheight = []; %Height of delta peak
    nearpeakdelta = zeros(2.*win.*deltaLFP.samplingRate+1,nummagbins);
    nearpeakgammaSPK = zeros(2.*win.*deltaLFP.samplingRate+1,nummagbins);
    DELTAmagbins = linspace(min(DELTAPeakheight),min([max(DELTAPeakheight)-1,4]), nummagbins);
    for pp = 1:length(sampleDELTA)
        ss = sampleDELTA(pp);
        switch NOSPIKES
            case false
                %Spikes around the delta
                nearpeakspikes = allspikes >= DELTApeaks(ss)-win & allspikes <= DELTApeaks(ss)+win;
                relspktime = [relspktime; allspikes(nearpeakspikes)-DELTApeaks(ss)];
                spkpeakheight = [spkpeakheight; DELTAPeakheight(ss).*ones(sum(nearpeakspikes),1)];
            case true
                %Gamma around the delta
                [~,groupidx] = min(abs(DELTAPeakheight(ss)-DELTAmagbins)); %Which row (magnitude) is this peak in?
                nearpeakgammaSPK(:,groupidx) =nansum([nearpeakgammaSPK(:,groupidx), ...
                    gammaLFP.normamp(DELTApeakIDX(ss)+[-1.*gammaLFP.samplingRate:gammaLFP.samplingRate])],2);
        end
            %Delta around the delta
            [~,groupidx] = min(abs(DELTAPeakheight(ss)-DELTAmagbins)); %Which row (magnitude) is this peak in?
            nearpeakdelta(:,groupidx) =nansum([nearpeakdelta(:,groupidx), ...
                deltaLFP.normamp(DELTApeakIDX(ss)+[-1.*deltaLFP.samplingRate:deltaLFP.samplingRate])],2);
    end
    peakmagdist = hist(DELTAPeakheight(sampleDELTA),DELTAmagbins);
    deltapower_byDELTAmag = bsxfun(@(A,B) A./B,nearpeakdelta,peakmagdist);
    switch NOSPIKES
        case false
            %Mean Spike rate around delta peaks
            spikehistmat = hist3([relspktime,spkpeakheight],{timebins,DELTAmagbins});
            ratemat_byDELTAmag = bsxfun(@(A,B) A./B./mean(diff(timebins))./numcells,spikehistmat,peakmagdist);
        case true
            %Mean High gamma around delta peaks
            ratemat_byDELTAmag = bsxfun(@(A,B) A./B,nearpeakgammaSPK,peakmagdist);
            lfptimebins = -win:1./deltaLFP.samplingRate:win;
            ratemat_byDELTAmag = interp1(lfptimebins,ratemat_byDELTAmag,timebins);
    end

    
    %GAMMA
    display('GAMMA...')
    %Pick the delta dips to average spiking around
    if length(GAMMAdips)>numchecks
        sampleGAMMA = randsample(length(GAMMAdips),numchecks); %Should really try to get even# in each bin...
    else
        sampleGAMMA = 1:length(GAMMAdips);
    end
    %Remove sample deltas too close to begining/end of recording
        sampleGAMMA(GAMMAdips(sampleGAMMA)<=win) = [];
        sampleGAMMA(GAMMAdips(sampleGAMMA)>=(gammaLFP.timestamps(end)-win)) = [];
    %Loop Through gamma dips and calculate spike rate (or high gamma)
    %around them, as a function of their peak height (= dip depth)
    relspktime = [];    %Spike time relative to delta peaks
    spkpeakheight = []; %Height of delta peak
    GAMMAmagbins = linspace(min(GAMMAdipdepth),min([max(GAMMAdipdepth)-0.5,2.5]),nummagbins);
    neardipgamma = zeros(2.*win.*gammaLFP.samplingRate+1,nummagbins);
    for pp = 1:length(sampleGAMMA)
        ss = sampleGAMMA(pp);
        switch NOSPIKES
            case false
        %Spikes around the gamma
        nearpeakspikes = allspikes >= GAMMAdips(ss)-win & allspikes <= GAMMAdips(ss)+win;
        relspktime = [relspktime; allspikes(nearpeakspikes)-GAMMAdips(ss)];
        spkpeakheight = [spkpeakheight; GAMMAdipdepth(ss).*ones(sum(nearpeakspikes),1)];
        end
         %Gamma around the gamma
        [~,groupidx] = min(abs(GAMMAdipdepth(ss)-GAMMAmagbins)); %Which gamma power bin
        neardipgamma(:,groupidx) =nansum([neardipgamma(:,groupidx), ... %nans for non-NREM from window smoothing
            gammaLFP.normamp(GAMMAdipIDX(ss)+[-1.*gammaLFP.samplingRate:gammaLFP.samplingRate])],2);
    end
    peakmagdist = hist(GAMMAdipdepth(sampleGAMMA),GAMMAmagbins);
    gammapower_byGAMMAmag = bsxfun(@(A,B) A./B,neardipgamma,peakmagdist);
    switch NOSPIKES
        case false
            %Mean Spike rate around gamma dips
    spikehistmat = hist3([relspktime,spkpeakheight],{timebins,GAMMAmagbins});
    ratemat_byGAMMAmag = bsxfun(@(A,B) A./B./mean(diff(timebins))./numcells,spikehistmat,peakmagdist);
        case true
            %Mean High gamma around gamma dips
        ratemat_byGAMMAmag = interp1(lfptimebins,gammapower_byGAMMAmag,timebins);
    end
    
    %% Find Thresholds as values where spike rate falls below threshold
    %Set the rate threshold to halfway between that at the delta peak and
    %the average around DELTA
    ratemat_byDELTAmag = imgaussfilt(ratemat_byDELTAmag,2);
    ratemat_byGAMMAmag = imgaussfilt(ratemat_byGAMMAmag,2);
    meanratearoundDELTA = nanmean(ratemat_byDELTAmag(:));
    minrateatDELTApeak = min(ratemat_byDELTAmag(round(end/2),:));
    DELTAraterange = meanratearoundDELTA-minrateatDELTApeak;
    
    ratethresh_Hz = minrateatDELTApeak+(ratethresh.*DELTAraterange);

    DELTAbox=bwmorph(ratemat_byDELTAmag<ratethresh_Hz,'close');
    DELTAbox=bwmorph(DELTAbox,'open');
    if sum(DELTAbox(:))== 0
        display('No DOWN around delta peak.... perhaps adjust rate threshold or pick another channel?')
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
        display('This may also indicate nan bug in delta rate threshold... DL ')
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
    %if SHOWFIG
        
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
    %end
%% Figure for lab meeting - illustrating threshold procedure
% xwin = [-0.5 0.5];
% exampledelta = [7,13,20];
% figure
% for ee = 1:length(exampledelta)
% subplot(6,3,ee+3)
%     plot(timebins,normrate_DELTA(:,exampledelta(ee)))
%     hold on
%     plot(xwin,ratethresh.*[1 1],'r--')
%     xlim(xwin);ylim([0 1.5])
%     ylabel('Mean Pop Rate (mean^-^1)');xlabel('t - aligned to delta peak (s)')
% subplot(6,3,ee)
%     plot(linspace(timebins(1),timebins(end),length(deltapower_byDELTAmag(:,exampledelta(ee)))),...
%         deltapower_byDELTAmag(:,exampledelta(ee)),'k','linewidth',2)
%     title(['Delta Peak: ',num2str(DELTAmagbins(exampledelta(ee)))])
%     xlim(xwin);ylim([-2 3.5])
%     ylabel('Delta-Filtered LFP');
% end
% NiceSave('ThresholdIllustration',basePath,baseName)
        
end


