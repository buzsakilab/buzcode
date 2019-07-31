function  [CellClass] = bz_CellClassification (basePath, varargin)
% Loads cell spike waveforms from the local folder, characterizes them and
% separates them into E vs I cells.  Manual verification based on clickable
% gui from Shige.
% 
%INPUTS
% baseName - basename of files in the local folder (default: pwd)
% 'knownE' - UIDs of known E cells, ie from synaptic interactions
% 'knownI' - UIDs of known I cells, ie from synaptic interactions
% 'ignorecells' - UIDs of cells to ignore in the clustering
% 'keepKnown' - keep the knownE/knownI labels (default: true) 
% 'saveMat'- true/false, save basePath/baseName.CellClass.cellinfo.mat
%            (default:true)
% 'saveFig'- true/false, save a DetectionFigure for posterity/QC 
%            (default:true)
% 'showFig'  show the figure without saving (default: false)
% 'forceReload'     -logical (default=false) to force reclassifying even if
%                    the CellClass.cellinfo.mat already exists
% 'noPrompts'          -logical (default) to supress any user prompts
%
%OUTPUTS
%   CellClass   buzcode structure saved to
%               basePath/baseName.CellClass.cellinfo.mat
%       .UID    -UID for each of the cells, matching spikes.cellinfo.mat
%       .pE 	-index vector, true for putative excitatory (RS) cells
%       .pI     -index vector, true for putative inhibitory (NS) cells
%       .label 	-labels for each cell 'pE' or 'pI'
%       .detectionparms.Waveforms -mean waveforms of each cell at the max channel
%       .detectionparms.TroughPeakMs
%       .detectionparms.SpikeWidthMs
%       .detectionparms.PyrBoundary - x,y of manually drawn line of boundary
%
%
% Mixture of functions from Shigeyoshi Fujisawa (WaveShapeClassification),
% Adrien Peyrache (wavelet-based determination of spike width) and Eran Stark (wfeatures, spikestats).
%
% Brendon Watson 2014
% Modified to buzcode format DLevenstein 2017
%% input Parsing

p = inputParser;
addParameter(p,'knownE',[],@isvector);
addParameter(p,'knownI',[],@isvector);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'showFig',false,@islogical);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'keepKnown',true,@islogical);
addParameter(p,'ignorecells',[],@isvector);

parse(p,varargin{:})

knownE = p.Results.knownE;
knownI = p.Results.knownI;
SAVEMAT = p.Results.saveMat;
SAVEFIG = p.Results.saveFig;
SHOWFIG = p.Results.showFig;
FORCERELOAD = p.Results.forceReload;
noPrompts = p.Results.noPrompts;
keepKnown = p.Results.keepKnown;
ignorecells = p.Results.ignorecells;
%%
Par = bz_getSessionInfo(basePath, 'noPrompts', noPrompts);
OneMs = round(Par.rates.wideband/1000);

baseName = bz_BasenameFromBasepath(basePath);
figfolder = fullfile(basePath,'DetectionFigures');
savefile = fullfile(basePath,[baseName,'.CellClass.cellinfo.mat']);
spikesfile = fullfile(basePath,[baseName,'.spikes.cellinfo.mat']);

if exist(savefile,'file') && ~FORCERELOAD
    display(['Cells already Classified, loading ',baseName,'.CellClass.cellinfo.mat'])
    load(savefile) %replace this with a bz_LoadCellinfo... function
    return
end
%% gather waves
if ~exist(spikesfile,'file')
    display(['spikes.cellinfo.mat does not yet exist,',...
        'saving one to insure cell UIDs are consistent across cellinfo files.'])
end
spikes = bz_GetSpikes('basepath',basePath,'saveMat',true);

MaxWaves = [];
for a = 1:length(spikes.rawWaveform)%in case stored wrong
    if iscolumn(spikes.rawWaveform{1})
        MaxWaves = cat(2,spikes.rawWaveform{:});
    elseif isrow(spikes.rawWaveform{1})
        MaxWaves = cat(1,spikes.rawWaveform{:})';
    else
        error('Something is wrong with your waveforms')
    end
end


%% get trough-peak delay times
%AllWaves(:,:,1) = [];
for a = 1:size(MaxWaves,2)
    thiswave = MaxWaves(:,a);
    [minval,minpos] = min(thiswave);
    minpos = minpos(1);
    [maxval,maxpos] = max(thiswave);
        [dummy,maxpos] = max(thiswave(minpos+1:end));
        if isempty(maxpos)
            warning('Your Waveform may be erroneous')
            maxpos = 1
        end
        maxpos=maxpos(1);
        maxpos = maxpos+minpos;
        tp(a) = maxpos-minpos; %In number of samples
end

%% get spike width by taking inverse of max frequency in spectrum (based on Adrien's use of Eran's getWavelet)
for a = 1:size(MaxWaves,2)
    w = MaxWaves(:,a);
    w = [w(1)*ones(1000,1);w;w(end)*ones(1000,1)];
    [wave f t] = getWavelet(w,20000,500,3000,128);
    %We consider only the central portion of the wavelet because we
    %haven't filtered it before hand (e.g. with a Hanning window)
    wave = wave(:,int16(length(t)/4):3*int16(length(t)/4));
    %Where is the max frequency?
    [maxPow ix] = max(wave);
    [dumy mix] = max(maxPow);
    ix = ix(mix);
    spkW(a) = 1000/f(ix);
end

%% Generate separatrix for cells 
x = tp'/OneMs;%trough to peak in ms
y = spkW';%width in ms of wavelet representing largest feature of spike complex... ie the full trough including to the tip of the peak

xx = [0 0.8];
yy = [2.4 0.4];
m = diff( yy ) / diff( xx );
b = yy( 1 ) - m * xx( 1 );  % y = ax+b
RS = y>= m*x+b;
INT = ~RS;

%% Convert knownE and knownI UIDs in to indices
knownEidx = ismember(spikes.UID,knownE);
knownIidx = ismember(spikes.UID,knownI);
ignoreidx = ismember(spikes.UID,ignorecells);

%% Plot for manual selection of boundary, with display of separatrix as a guide.

if all(knownEidx | knownIidx | ignoreidx) && keepKnown
    PyrBoundary = [nan nan];
    ELike = false(size(spikes.UID));
else
    h = figure('position',[674   456   560   420]);
    fprintf('\nDiscriminate pyr and int (select Pyramidal)');
    xlabel('Trough-To-Peak Time (ms)')
    ylabel('Wave width (via inverse frequency) (ms)')
    figure('position',[674   961   561   109]);
    title({'Discriminate pyr and int (select Pyramidal)','left click to draw boundary', 'center click/ENTER to complete)'});
    figure(h)
    [ELike,PyrBoundary] = ClusterPointsBoundaryOutBW([x y],knownEidx,knownIidx,m,b);
    ELike = ELike';
end

if keepKnown
    ELike(knownEidx) = 1;
    ELike(knownIidx) = 0;
end
%% Mean waveforms output
CellClass.UID = spikes.UID;
CellClass.pE = ELike & ~ignoreidx;
CellClass.pI = ~ELike & ~ignoreidx;
CellClass.label = cell(size(CellClass.UID));
CellClass.label(CellClass.pE) = {'pE'};
CellClass.label(CellClass.pI) = {'pI'};
CellClass.celltypes = {'pE','pI'}; %should this include pI if there are no pI cells?
CellClass.detectionparms.TroughPeakMs = x';
CellClass.detectionparms.SpikeWidthMs = y';
CellClass.detectionparms.PyrBoundary = PyrBoundary;
CellClass.detectionparms.Waveforms = MaxWaves;

if SAVEMAT
    save(savefile,'CellClass')
end

%%
if SAVEFIG || SHOWFIG
    figure
    subplot(2,2,1)
        plot(CellClass.detectionparms.TroughPeakMs(CellClass.pE),...
            CellClass.detectionparms.SpikeWidthMs(CellClass.pE),'.','color',[0 0.6 0])
        hold on
        plot(CellClass.detectionparms.TroughPeakMs(CellClass.pI),...
            CellClass.detectionparms.SpikeWidthMs(CellClass.pI),'.','color',[0.6 0 0])
        plot(CellClass.detectionparms.TroughPeakMs(ignoreidx),...
            CellClass.detectionparms.SpikeWidthMs(ignoreidx),'x','color',[0.5 0.5 0.5])
        axis tight
        plot(CellClass.detectionparms.PyrBoundary(:,1),...
            CellClass.detectionparms.PyrBoundary(:,2))
        xlim([0 max([x+0.1;2])])
        ylim([0 max([y+0.1;2])])
        xb = get(gca,'XLim');
        yb = get(gca,'YLim');
        plot(xb,[m*xb(1)+b m*xb(2)+b])
        xlabel('Trough to Peak Time (ms)')
        ylabel('Spike Width (ms)')
        title([baseName,': Cell Classification'])
        
    subplot(2,2,2)
        if any(CellClass.pE)
            plot([1:size(MaxWaves,1)]./OneMs,MaxWaves(:,CellClass.pE),'color',[0 0.6 0])
        end
        hold on
        if any(CellClass.pI)
            plot([1:size(MaxWaves,1)]./OneMs,MaxWaves(:,CellClass.pI),'color',[0.6 0 0])
        end
        if any(ignoreidx)
            plot([1:size(MaxWaves,1)]./OneMs,MaxWaves(:,ignoreidx),'color',[0.5 0.5 0.5])
        end
        axis tight
        xlabel('t (ms)')
        if SAVEFIG
        NiceSave('CellClassification',figfolder,baseName)
        end
end
