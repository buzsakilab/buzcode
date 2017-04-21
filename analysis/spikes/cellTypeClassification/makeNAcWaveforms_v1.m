% this is a modified version of Brendon's cell type classification code... both should be merged into a single file
function makeNAcWaveforms_v1(basedir, cellshanknums, shanks)

% function makeNAcWaveforms_v1(basedir, cellshanknums, shanks)
%
% A modified version of Brendon's code. - Luke
%
% Loads cell spike waveforms from the local folder, characterizes them and
% separates them into E vs I cells.  Manual verification based on clickable
% gui from Shige.
% 
% INPUTS
% filebasename - basename of files in the local folder
% cellshanknums - per-shank number of each cell to load, ie cell 2 on shank 4 has a
%                 cellshanknum of 2
% shanks - spikegroup numbers of each cell, these corresond to .clu.x
%                 numbers
% ESynapseCells - pre-labeled E cells, ie from synaptic interactions
% ISynapseCells - pre-labeled I cells, ie from synaptic interactions
%
% OUTPUTS
% CellClassOutput - array with multiple columns:
%        - Column 1: cell numbers in order
%        - Column 2: trough-peak lag in ms
%        - Column 3: wavelet-based spike width
%        - Column 4: neurons classified as having I-like waveforms
%        - Column 5: neurons classified as having E-like waveforms
% PyrBoundary - x,y of manually drawn line of boundary
% Waveforms - the saved waveforms per cell
%
% Mixture of functions from Shigeyoshi Fujisawa (WaveShapeClassification),
% Adrien Peyrache (wavelet-based determination of spike width) and Eran Stark (wfeatures, spikestats).
%
% Brendon Watson 2014

cd(basedir);
cd('DAT');

[~, filebasename, ~] = fileparts(pwd);
Par = LoadPar([filebasename '.xml']);
cellnums = 1:length(shanks);

OneMs = round(Par.SampleRate/1000);
MaxWaves = [];
allshanks = unique(shanks);

%% gather waves
for a = 1:length(allshanks);
    thisshank = allshanks(a);
    AllWaves{thisshank} = [];%separated AllWaves in case each shank has diff number of sites
%     theseChannels = Par.SpkGrps(thisshank).Channels;
    numChannels = length(Par.SpkGrps(thisshank).Channels);
    nSamples = Par.SpkGrps(a).nSamples;
    spkname = fullfile([basedir '/DAT'],[filebasename '.spk.' num2str(thisshank)]);
    
    cellsthisshank = cellnums(shanks==thisshank);
    intrashankclunums = cellshanknums(shanks == thisshank);
    cluname = fullfile([basedir '/DAT'],[filebasename '.clu.' num2str(thisshank)]);
    clu = LoadClu(cluname);
    for b = 1:length(cellsthisshank)
        spikesthiscell = find(clu == intrashankclunums(b));        
        Waveforms = LoadSpikeWaveforms_BW(spkname,numChannels,nSamples,spikesthiscell);%load one cell at a time
        meanwaves = squeeze(mean(Waveforms,1));
        AllWaves{thisshank}(:,:,end+1) = meanwaves;
        [dummy,maxwaveidx] = max(abs(max(meanwaves,[],2)-min(meanwaves,[],2)));
        MaxWaves(:,end+1) = meanwaves(maxwaveidx,:);
    end
    disp(['Shank ',num2str(a),' (Orig#:' num2str(thisshank) ') Done'])
    AllWaves{thisshank}(:,:,1) = [];
end

%% get trough-peak delay times
AllWaves(:,:,1) = [];
for a = 1:length(cellnums)
    thiswave = MaxWaves(:,a);
    [minval,minpos] = min(thiswave);
    minpos = minpos(1);
    [maxval,maxpos] = max(thiswave);
        [dummy,maxpos] = max(thiswave(minpos+1:end));
        maxpos=maxpos(1);
        maxpos = maxpos+minpos;
        tp(a) = maxpos-minpos;
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



%% Mean waveforms output
TroughPeakMs = x;
SpikeWidthMs = y;
NAcWaveforms_v1 = v2struct(AllWaves,MaxWaves,cellnums,cellshanknums,shanks,filebasename,TroughPeakMs,SpikeWidthMs);
NAcWaveforms_v1.basedir = basedir;

save([basedir '/DAT/Analysis/NAcWaveforms_v1.mat'], 'NAcWaveforms_v1');

