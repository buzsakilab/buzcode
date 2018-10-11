function  [CellClass] = bz_CellClassification (basePath, varargin)
% Loads cell spike waveforms from the local folder, characterizes them and
% separates them into E vs I cells.  Manual verification based on clickable
% gui from Shige.
% 
%INPUTS
% baseName - basename of files in the local folder (default: pwd)
% 'knownE' - UIDs of known E cells, ie from synaptic interactions
% 'knownI' - UIDs of known I cells, ie from synaptic interactions
% 'saveMat'- true/false, save basePath/baseName.CellClass.cellinfo.mat
%            (default:true)
% 'saveFig'- true/false, save a DetectionFigure for posterity/QC 
%            (default:true)
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
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'noPrompts',false,@islogical);

parse(p,varargin{:})

knownE = p.Results.knownE;
knownI = p.Results.knownI;
SAVEMAT = p.Results.saveMat;
SAVEFIG = p.Results.saveFig;
FORCERELOAD = p.Results.forceReload;
noPrompts = p.Results.noPrompts;
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
MaxWaves = cat(1,spikes.rawWaveform{:})';

%% Previous waveform loading code
% for a = 1:length(allshanks);
%     thisshank = allshanks(a);
%     AllWaves{thisshank} = [];%separated AllWaves in case each shank has diff number of sites
% %     theseChannels = Par.SpkGrps(thisshank).Channels;
%     numChannels = length(Par.SpkGrps(thisshank).Channels);
%     nSamples = Par.SpkGrps(a).nSamples;
%     spkname = fullfile(fpath,[filebasename '.spk.' num2str(thisshank)]);
%     
%     cellsthisshank = cellnums(shanks==thisshank);
%     intrashankclunums = cellshanknums(shanks == thisshank);
%     cluname = fullfile(fpath,[filebasename '.clu.' num2str(thisshank)]);
%     clu = LoadClu(cluname);
%     for b = 1:length(cellsthisshank)
%         spikesthiscell = find(clu == intrashankclunums(b));        
%         Waveforms = LoadSpikeWaveforms(spkname,numChannels,nSamples,spikesthiscell);%load one cell at a time
%         meanwaves = squeeze(mean(Waveforms,1));
%         AllWaves{thisshank}(:,:,end+1) = meanwaves;
%         [dummy,maxwaveidx] = max(abs(max(meanwaves,[],2)-min(meanwaves,[],2)));
%         MaxWaves(:,end+1) = meanwaves(maxwaveidx,:);
%     end
%     disp(['Shank ',num2str(a),' (Orig#:' num2str(thisshank) ') Done'])
%     AllWaves{thisshank}(:,:,1) = [];
% end

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

%% Plot for manual selection of boundary, with display of separatrix as a guide.
h = figure;
title({'Discriminate pyr and int (select Pyramidal)','left click to draw boundary', 'center click/ENTER to complete)'});
fprintf('\nDiscriminate pyr and int (select Pyramidal)');
xlabel('Trough-To-Peak Time (ms)')
ylabel('Wave width (via inverse frequency) (ms)')
[ELike,PyrBoundary] = ClusterPointsBoundaryOutBW([x y],knownE,knownI,m,b);

%% Mean waveforms output
CellClass.UID = spikes.UID;
CellClass.pE = ELike';
CellClass.pI = ~ELike';
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
if SAVEFIG
    figure
    subplot(2,2,1)
        plot(CellClass.detectionparms.TroughPeakMs(CellClass.pE),...
            CellClass.detectionparms.SpikeWidthMs(CellClass.pE),'k.')
        hold on
        plot(CellClass.detectionparms.TroughPeakMs(CellClass.pI),...
            CellClass.detectionparms.SpikeWidthMs(CellClass.pI),'r.')
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
        plot([1:size(MaxWaves,1)]./OneMs,MaxWaves(:,CellClass.pE),'color',[0 0.6 0])
        hold on
        plot([1:size(MaxWaves,1)]./OneMs,MaxWaves(:,CellClass.pI),'color',[0.6 0 0])
        axis tight
        xlabel('t (ms)')
        
        NiceSave('CellClassification',figfolder,baseName)
end
