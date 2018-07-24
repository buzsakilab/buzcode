
function [spikes] = bz_LoadPhy(varargin)
% Load kilosort/phy clusters
%
% USAGE
%
%   [spikes] = bz_loadPhy(varargin)
%
% INPUTS
%
% basepath        -path to recording (where .dat are)
% kilosort_path   -path to kilosort* folder after phy curing (by default
%   .\kilosort_\)
% getWaveforms    -logical (default=true) to load mean of raw waveform data
% saveMat         -logical (default=false) to save in buzcode format
% UID             -vector subset of UID's to load 
% fs              -scalar (default is taken from sessionInfo, otherwise
%   30000). Sampling frequency. 
% nChannels       -scalar (default is taken from sessionInfo, otherwise
%   32). Total number of channels.
% forceReload    -logical (default=false) to force loading from
%                     res/clu/spk files
%
% OUTPUTS
%
% spikes - cellinfo struct with the following fields
%   .sessionName    -name of recording file
%   .UID            -unique identifier for each neuron in a recording
%   .times          -cell array of timestamps (seconds) for each neuron
%   .spindices      -sorted vector of [spiketime UID], useful for 
%                           input to some functions and plotting rasters
%   .shankID        -shank ID that each neuron was recorded on
%   .amplitude:     amplitudes for cluster N
%   .maxWaveformCh  -channel # with largest amplitude spike for each neuron
%   .rawWaveform    -average waveform on maxWaveformCh (from raw .dat)
%
%   Manu Valero 2018

% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'kilosort_path',ls('Kilosort*'),@isstr); % probably this line only works in windows
addParameter(p,'getWaveforms',true,@islogical)
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'UID',[],@isvector);
addParameter(p,'fs',30000,@isnumeric);
addParameter(p,'nChannels',32,@isnumeric);
addParameter(p,'forceReload',false,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
kilosort_path = p.Results.kilosort_path;
getWave = p.Results.getWaveforms;
saveMat = p.Results.saveMat;
UID = p.Results.UID;
fs = p.Results.fs; % it will be overwritten if bz_getSessionInfo
nChannels = p.Results.nChannels; % it will be overwritten if bz_getSessionInfo
forceReload = p.Results.forceReload;

try [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', false);
    fs = sessionInfo.rates.wideband;
    nChannels = sessionInfo.nChannels;
catch
    warning('SessionInfo file not found.');
end

if exist([basepath filesep sessionInfo.FileName '.spikes.cellinfo.mat'],'file') && forceReload == false
    disp('loading spikes from cellinfo file..')
    load([basepath filesep sessionInfo.FileName '.spikes.cellinfo.mat'])
    
else
    spike_cluster_index = readNPY(fullfile(kilosort_path, 'spike_clusters.npy'));
    spike_times = readNPY(fullfile(kilosort_path, 'spike_times.npy'));
    % spike_amplitudes = readNPY(fullfile(kilosort_path, 'amplitudes.npy'));
    % spike_clusters = unique(spike_cluster_index);
    cluster_group = tdfread(fullfile(kilosort_path,'cluster_group.tsv'));
    shanks = readNPY(fullfile(kilosort_path, 'shanks.npy')); % done

    % parameters for extract raw waveforms
    if getWave
        gwfparams.dataDir = strcat(pwd,'\',kilosort_path);    % KiloSort/Phy output folder
        basepath = cd;
        [~,basename] = fileparts(basepath);
        gwfparams.fileName = [basepath,'\',basename,'.dat']; % AP band file from spikeGLX specifically
        gwfparams.dataType = 'int16';            % Data type of 000.dat file (this should be BP filtered)
        gwfparams.nCh = nChannels;               % Number of channels that were streamed to disk in .dat file
        gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
        gwfparams.nWf = 2000;                      % Number of waveforms per unit to pull out
    end

    spikes = [];
    spikes.sessionName = sessionInfo.FileName;
    jj = 1;
    for ii = 1:length(cluster_group.group)
        if strcmp(cluster_group.group(ii,:),'good ')
            ids = find(spike_cluster_index == cluster_group.cluster_id(ii)); % cluster id
            spikes.UID(jj) = cluster_group.cluster_id(ii);
            spikes.times{jj} = double(spike_times(ids))/fs; % cluster time
            spikes.ts{jj} = double(spike_times(ids)); % cluster time
            [~,cluster_id] = find(cluster_group.cluster_id == spikes.UID(jj));
            spikes.shankID(jj) = shanks(cluster_id);
            % spikes.amplitudes{jj} = double(spike_amplitudes(ids));

            if getWave
                gwfparams.spikeTimes = spike_times(ids);
                gwfparams.spikeClusters = spike_cluster_index(ids);
                wf = getWaveForms(gwfparams);
                wf.waveFormsMean = squeeze(wf.waveFormsMean);
                [~, spikes.maxWaveformCh(jj)] = max(wf.waveFormsMean(:,41));
                spikes.rawWaveform{jj} = wf.waveFormsMean(spikes.maxWaveformCh(jj),:);
            end

            jj = jj + 1;
            if jj > 2
                fprintf(repmat('\b', 1, 9));
            end
            fprintf(' %2.1f%%...',ii/length(cluster_group.group)*100);
        end
    end
    fprintf('\n');

    % saveMat (only saving if no exclusions)
    if saveMat
        save([basepath filesep sessionInfo.FileName '.spikes.cellinfo.mat'],'spikes');
    end
end

% filter by UID input
if ~isempty(UID)
      [toRemove] = ~ismember(spikes.UID,UID);
    spikes.UID(toRemove) = [];
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.times{r} = [];
         % spikes.region{r} = [];
        end
    end
    spikes.times = removeEmptyCells(spikes.times);
    %spikes.region = removeEmptyCells(spikes.region);
    % spikes.cluID(toRemove) = [];
    spikes.shankID(toRemove) = [];
    
    if getWaveforms
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.rawWaveform{r} = [];
        end
    end
    spikes.rawWaveform = removeEmptyCells(spikes.rawWaveform);
    spikes.maxWaveformCh(toRemove) = [];
    end
end

% Generate spindices matrics
for cc = 1:length(spikes.UID)
    groups{cc}=spikes.UID(cc).*ones(size(spikes.times{cc}));
end
if ~isempty(spikes.UID)
    alltimes = cat(1,spikes.times{:}); groups = cat(1,groups{:}); %from cell to array
    [alltimes,sortidx] = sort(alltimes); groups = groups(sortidx); %sort both
    spikes.spindices = [alltimes groups];
end

% Compute spike measures
if ~isempty(spikes.UID)
    for ii = 1:size(spikes.UID,2)
        [~,tmp] = max(spikes.rawWaveform{ii}(1,41:end));
        p2pWidth = tmp/fs; % peak (negative) to peak (second positive) duration 
        [spikes.autocorr{ii}.xout,spikes.autocorr{ii}.r,spikes.autocorr{ii}.peakAutocorr] =...
            autocorr_spikes(spikes.ts{ii},fs,26,1);
    end 
end