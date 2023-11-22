function [ISILFPMap] = bz_ISILFPMap(basePath,varargin)
% Under Development: DLevenstein 2020
%% Parse the inputs
defaultstates.ALL = [-Inf Inf];

p = inputParser;
addParameter(p,'ints',defaultstates)
addParameter(p,'savechaninfo',false,@islogical)
addParameter(p,'cellclass','AllCells');
addParameter(p,'regions','sessionInfo');
addParameter(p,'groups','sessionInfo');
addParameter(p,'figfolder',false)
addParameter(p,'figname',[])
addParameter(p,'showfig',false,@islogical);
addParameter(p,'forceRedetect',false,@islogical);
addParameter(p,'dt',0.25);
addParameter(p,'winsize',1);
addParameter(p,'frange',[1 312]);
addParameter(p,'nfreqs',150);

parse(p,varargin{:})
ints = p.Results.ints;
cellclass = p.Results.cellclass;
SAVECHANINFO = p.Results.savechaninfo;
regions = p.Results.regions;
groups = p.Results.groups;
figfolder = p.Results.figfolder;
figname = p.Results.figname;
SHOWFIG = p.Results.showfig;
forceRedetect = p.Results.forceRedetect;
dt = p.Results.dt;
winsize = p.Results.winsize;
frange = p.Results.frange;
nfreqs = p.Results.nfreqs;

%% Load Header
%Initiate Paths
%reporoot = '/Users/dl2820/Project Repos/NeuronalHeterogeneity/';
% basePath = '/Users/dl2820/Dropbox/Research/Datasets/20140526_277um';
%basePath = '/Users/dl2820/Dropbox/Research/Datasets/Cicero_09102014';
% %basePath = pwd;
% %basePath = fullfile(reporoot,'Datasets/onProbox/AG_HPC/Achilles_11012013');
%figfolder = [reporoot,'AnalysisScripts/AnalysisFigs/DailyAnalysis'];
baseName = bz_BasenameFromBasepath(basePath);
chaninfofilename = fullfile(basePath,[baseName,'.ISILFPMap.channelinfo.mat']);

if exist(chaninfofilename,'file') && ~forceRedetect
    ISILFPMap = load(chaninfofilename);
    return
end

%Load Stuff
sessionInfo = bz_getSessionInfo(basePath);
spikes = bz_GetSpikes('basePath',basePath,'noPrompts',true);

switch cellclass
    case 'load'
        CellClass = bz_LoadCellinfo(basePath,'CellClass');
    case 'AllCells'
        CellClass.AllCells = true(size(spikes.times));
        CellClass.celltypes = {'AllCells'};
end
try
celltypes = CellClass.celltypes;
catch
    celltypes = unique(CellClass.label);
end
cellcolor = {'k','r'};

if strcmp(ints,'load');
        SleepState = bz_LoadStates(basePath,'SleepState');
        ints = SleepState.ints;
end
states = fieldnames(ints);
%statecolors = {[0 0 0],[0 0 1],[1 0 0]};



%%



%%
if strcmp(groups,'sessionInfo')
        groups = sessionInfo.spikeGroups.groups{:};
end
if strcmp(regions,'sessionInfo')
        try
            emptychans = cellfun(@isempty,sessionInfo.region);
            Regions = unique(sessionInfo.region(~emptychans));
            sessionInfo.region(emptychans) = {'empty'};
        catch
            display('No regions in sessionInfo, using all channels in spikegroups')
            Regions{1} = 'NA';
            inSGchans =[groups{:}];
            inSGchans = ismember(sessionInfo.channels,inSGchans);
            sessionInfo.region(inSGchans) = {'NA'};
            
            emptychans = cellfun(@isempty,sessionInfo.region);
            sessionInfo.region(emptychans) = {'empty'};
        end
end
%Regions = unique(spikes.region(~cellfun(@isempty,spikes.region)));

ISILFPMap.Regions = Regions;
%% Load the LFP
for rr = 1:length(Regions)
    display(['Region: ',Regions{rr}])
    inregionchanIDX = ismember(sessionInfo.region,Regions{rr});
    inregionchan = sessionInfo.channels(inregionchanIDX);
    if isempty(inregionchan)
        display('No LFP Channels')
        continue
    end
    try
        inregioncellIDX = ismember(spikes.region,Regions{rr}) | ismember(spikes.maxWaveformCh,inregionchan);
    catch
        display('No Regions in spikes.cellinfo.mat')
        inregioncellIDX = ismember(spikes.maxWaveformCh,inregionchan);
    end
    if sum(inregioncellIDX)==0
        display('No Cells')
        continue
    end
    downsamplefactor = 1;
    lfp = bz_GetLFP(inregionchan,... %note: may have to load separately here for RAM...
        'basepath',basePath,'noPrompts',true,'downsample',downsamplefactor);

%% Calculate Residuals
[specslope] = bz_PowerSpectrumSlope(lfp,winsize,dt,'spectype','fft',...
    'nfreqs',nfreqs,'showfig',false,'showprogress',false,'frange',frange);%,...
    %'saveMat',basePath,'saveName',['Chan',num2str(lfpchannel)],...
    %'saveFolder','WavPSS');
    

%%
ISILFPMap.(Regions{rr}).UIDs = spikes.UID(inregioncellIDX);
ISILFPMap.(Regions{rr}).chanID = inregionchan;
ISILFPMap.freqs = specslope.freqs;
ISILFPMap.baseName = baseName;
%%
for tt = 1:length(celltypes)
    popspikes.(Regions{rr}).(celltypes{tt}) = cat(1,spikes.times{CellClass.(celltypes{tt})&inregioncellIDX});
end

%%
clear fPower
%Do inintervals outside the loop, for speed
for ss = 1:length(states)
    fPower.(states{ss}).timestamps = specslope.timestamps;
    fPower.(states{ss}).idx = InIntervals(fPower.(states{ss}).timestamps,ints.(states{ss}));
    fPower.(states{ss}).timestamps = fPower.(states{ss}).timestamps(fPower.(states{ss}).idx);
end
%%

clear PopConditionalISIDist PopConditionalISIDist_phase PopConditionalISIDist_power
for ff = 1:length(specslope.freqs)
	bz_Counter(ff,length(specslope.freqs),'Freq')
    for cc = 1:length(inregionchan)
        fPower.data = specslope.resid(:,ff,cc);
        for tt = 1:length(celltypes)
            for ss = 1:length(states)
                %Nan Pad data to avoid having to use InIntervals
                fPower.(states{ss}).data = fPower.data(fPower.(states{ss}).idx);
                [ fPower.(states{ss}).data ] = NanPadJumps( fPower.(states{ss}).timestamps,...
                    fPower.(states{ss}).data,2*dt);
                %tic
                [PopConditionalISIDist] = bz_ConditionalISI(popspikes.(Regions{rr}).(celltypes{tt}),fPower.(states{ss}),...
                    'showfig',false,'ISIDist',false);%, 'ints',ints.(states{ss}));%,...
                %test(rr,ss,tt,ff,cc) = PopConditionalISIDist.MutInf;
                ISILFPMap.(Regions{rr}).(states{ss}).(celltypes{tt})(ff,cc) = PopConditionalISIDist.MutInf;

                [UnitsConditionalISIDist] = bz_ConditionalISI(spikes.times,fPower.(states{ss}),...
                    'showfig',false,'ISIDist',false);%,'ints',ints.(states{ss}));%,...

                ISILFPMap.(Regions{rr}).(states{ss}).AllUnits(ff,cc,:) = UnitsConditionalISIDist.MutInf;
                %toc
            end
        end
    end
    
end

%% Mean over all cells in region/cell type
for ff = 1:length(specslope.freqs)
    for cc = 1:length(inregionchan)
        for tt = 1:length(celltypes)
            for ss = 1:length(states)
                ISILFPMap.(Regions{rr}).(states{ss}).Mean.(celltypes{tt})(ff,cc) = ...
                    nanmean(ISILFPMap.(Regions{rr}).(states{ss}).AllUnits(ff,cc,CellClass.(celltypes{tt})&inregioncellIDX));
            end 
        end
    end
end



%%
%rr = 2;
[~,SGsIDX] = cellfun(@(SGAll) (ismember(SGAll,inregionchan)),...
                   groups,'UniformOutput',false); %Just cells in the right region
SGsThisRegion = cellfun(@(SGAll) SGAll(ismember(SGAll,inregionchan)),...
                   groups,'UniformOutput',false); %Just cells in the right region
SGLength = cellfun(@(SG) length(SG),SGsThisRegion);
SGnum = find(SGLength>0);
SGorder = [SGsIDX{:}]; %The spikegroup order of all cells 
SGorder = SGorder(SGorder~=0);
%SGsort = 
%Sort by group
%save MI map - make sure it has channel numbers and UID of cells in each
%region/population. and freqs. and spike group ID/sortings.
%Save figures in detection figures
ISILFPMap.(Regions{rr}).SGorder = SGorder;
ISILFPMap.(Regions{rr}).SGLength = SGLength;
ISILFPMap.(Regions{rr}).SGnum = SGnum;
%%
if SHOWFIG | figfolder
figure
    for ss = 1:length(states)
        for tt = 1:length(celltypes)
            if tt>2
               display(celltypes{tt})
               continue
            end
            subplot(2,3,ss+(tt-1).*3)
                imagesc(log2(ISILFPMap.freqs),[1 length(ISILFPMap.(Regions{rr}).chanID)],...
                    ISILFPMap.(Regions{rr}).(states{ss}).(celltypes{tt})(:,SGorder)')
                hold on
                plot(log2(ISILFPMap.freqs([1 end])),cumsum(SGLength)'*[1 1]+0.5,'w')
                
                set(gca,'YTickLabel',[]);
                if ss == 1
                    ylabel({(celltypes{tt}),'Spike Group'})
                    set(gca,'YTick',cumsum(SGLength(SGnum))-0.5.*SGLength(SGnum)+0.5)
                    set(gca,'YTickLabel',SGnum)
                end
                if tt == 1
                    title(states{ss})
                end
                if tt == length(celltypes)
                    xlabel('f (Hz)');
                    LogScale('x',2)
                else
                    set(gca,'XTickLabel',[]);
                end
        end
    end
    
if figfolder
    if ~isempty(figname)
        baseName = figname;
    end
    NiceSave(['ISIMod_',Regions{rr}],figfolder,baseName)
end
    

end
end
%%
if SAVECHANINFO
    save(chaninfofilename,'ISILFPMap')
end

end