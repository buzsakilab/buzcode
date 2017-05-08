function spikes = bz_GetSpikes(varargin)
% bz_getSpikes - Get spike timestamps.
%
% USAGE
%
%    spikes = bz_getSpikes(varargin)
% 
% INPUTS
%
%    spikeGroups     -vector subset of shank IDs to load
%    region          -string region ID to load neurons from specific region (requires metadata file)
%    UID             -vector subset of UID's to load 
%    basepath        -path to recording (where .dat/.clu/etc files are)
%    getWaveforms    -logical (default=true) to load mean of raw waveform data
%    forceReload     -logical (default=false) to force loading from
%                     res/clu/spk files
%    saveMat         -logical (default=false) to save in buzcode format
%    
% OUTPUTS
%
%    spikes - cellinfo struct with the following fields
%          .sessionName    -name of recording file
%          .UID            -unique identifier for each neuron in a recording
%          .times          -cell array of timestamps (seconds) for each neuron
%          .region         -region ID for each neuron (especially important large scale, high density probes)
%          .shankID        -shank ID that each neuron was recorded on
%          .maxWaveformCh  -channel # with largest amplitude spike for each neuron
%          .rawWaveform    -average waveform on maxWaveformCh (from raw .dat)
%          .cluID          -cluster ID, NOT UNIQUE ACROSS SHANKS
%           
% NOTES
%
% This function can be used in several ways to load spiking data.
% Specifically, it loads spiketimes for individual neurons and other
% metadata that describes each neuron.  Spiketimes can be loaded using the
% UID(1-N), the shank the neuron was on, or the region it was recorded in.
% The default behavior is to load all spikes in a recording. The .shankID
% and .cluID fields can be used to reconstruct the 'units' variable often
% used in FMAToolbox.
% units = [spikes.shankID spikes.cluID];
% 
% 
% first usage recommendation:
% 
%   spikes = bz_getSpikes('saveMat',true); Loads and saves all spiking data
%                                          into buzcode format .cellinfo. struct
% other examples:
%
%   spikes = bz_getSpikes('spikeGroups',1:5); first five shanks
%
%   spikes = bz_getSpikes('region','CA1'); cells tagged as recorded in CA1
%
%   spikes = bz_getSpikes('UID',[1:20]); first twenty neurons
%
%
% written by David Tingley, 2017

% TODO
% - integrate session metadata in place of xml-LoadParameters
% - get 'region' input working with session metadata
p = inputParser;
addParameter(p,'spikeGroups',[],@isvector);
addParameter(p,'region','',@isstr); % won't work without metadata 
addParameter(p,'UID',[],@isvector);
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'getWaveforms',true,@islogical)
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'saveMat',false,@islogical);

parse(p,varargin{:})

spikeGroups = p.Results.spikeGroups;
region = p.Results.region;
UID = p.Results.UID;
basepath = p.Results.basepath;
getWaveforms = p.Results.getWaveforms;
forceReload = p.Results.forceReload;
saveMat = p.Results.saveMat;

% get meta info about recording, for now we'll use the xml.
meta = LoadParameters(basepath);
samplingRate = meta.rates.wideband;
nChannels = meta.nChannels;


% to use session metadata, ADD IT HERE AS 'meta' and not below.

% samplingRate = meta.rates.wdieband;
% nChannels = meta.nChannels;


%% if the cellinfo file exist and we don't want to re-load files
if exist([basepath filesep meta.FileName '.spikes.cellinfo.mat']) & forceReload == false
    disp('loading spikes from cellinfo file..')
    load([basepath filesep meta.FileName '.spikes.cellinfo.mat'])
else % do the below then filter by inputs...
    
disp('loading spikes from clu/res/spk files..')
% find res/clu/fet/spk files here
cluFiles = dir([basepath filesep '*.clu*']);  % sort to be in correct order?
resFiles = dir([basepath filesep '*.res*']);
spkFiles = dir([basepath filesep '*.spk*']);

% ensures we load in sequential order (forces compatibility with FMAT
% ordering)
for i = 1:length(cluFiles)
    temp = strsplit(cluFiles(i).name,'.');
    shanks(i) = str2num(temp{length(temp)});
end
[shanks ind] = sort(shanks);
cluFiles = cluFiles(ind);
resFiles = resFiles(ind); 
spkFiles = spkFiles(ind);

% check if there are matching #'s of files
if length(cluFiles) ~= length(resFiles) & length(cluFiles) ~= length(spkFiles)
    error('found an incorrect number of res/clu/spk files...')
end

% use the .clu files to get spike ID's and generate UID and spikeGroup
% use the .res files to get spike times
count = 1;

for i=1:length(cluFiles) 
    disp(['working on ' spkFiles(i).name])
    
    temp = strsplit(cluFiles(i).name,'.');
    shankID = str2num(temp{length(temp)});
    clu = load(fullfile(basepath,cluFiles(i).name));
    clu = clu(2:end); % toss the first sample to match res/spk files
    res = load(fullfile(basepath,resFiles(i).name));
    nSamples = meta.spikeGroups.nSamples(shankID);
    spkGrpChans = meta.spikeGroups.groups{shankID}; % we'll eventually want to replace these two lines
    
    if getWaveforms
        % load waveforms
        chansPerSpikeGrp = length(meta.spikeGroups.groups{shankID});
        fid = fopen(fullfile(basepath,spkFiles(i).name),'r');
        wav = fread(fid,[1 inf],'int16');
        wav = reshape(wav,chansPerSpikeGrp,nSamples,[]);
        wav = permute(wav,[3 1 2]);
    end
    
    cells  = unique(clu);
    % remove MUA and NOISE clusters...
    cells(cells==0) = [];
    cells(cells==1) = [];  % consider adding MUA as another input argument...?
    
    for c = 1:length(cells)
       spikes.UID(count) = count; % this only works if all shanks are loaded... how do we optimize this?
       ind = find(clu == cells(c));
       spikes.times{count} = res(ind) ./ samplingRate;
       spikes.shankID(count) = shankID;
       spikes.cluID(count) = cells(c);
              
       if getWaveforms
           wvforms = squeeze(mean(wav(ind,:,:)))-mean(mean(mean(wav(ind,:,:)))); % mean subtract to account for slower (theta) trends
           for t = 1:size(wvforms,1)
              [a(t) b(t)] = max(abs(wvforms(t,:))); 
           end
           [aa bb] = max(a);
           spikes.rawWaveform{count} = wvforms(bb,:);
           spikes.maxWaveformCh(count) = spkGrpChans(bb);     
           clear a aa b bb
       end
       count = count + 1;
    end
end

spikes.sessionName = meta.FileName;

end


%% filter by spikeGroups input
if ~isempty(spikeGroups)
    [toRemove] = ~ismember(spikes.shankID,spikeGroups);
    spikes.UID(toRemove) = [];
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.times{r} = [];
        end
    end
    spikes.times = removeEmptyCells(spikes.times);
    spikes.cluID(toRemove) = [];
    spikes.shankID(toRemove) = [];
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.rawWaveform{r} = [];
        end
    end
    spikes.rawWaveform = removeEmptyCells(spikes.rawWaveform);
    spikes.maxWaveformCh(toRemove) = [];
end
%% filter by region input
if ~isempty(region)
    warning('region input is not working yet, get metadata loading working first!')
        [toRemove] = ~ismember(spikes.region,region);
    spikes.UID(toRemove) = [];
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.times{r} = [];
        end
    end
    spikes.times = removeEmptyCells(spikes.times);
    spikes.cluID(toRemove) = [];
    spikes.shankID(toRemove) = [];
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.rawWaveform{r} = [];
        end
    end
    spikes.rawWaveform = removeEmptyCells(spikes.rawWaveform);
    spikes.maxWaveformCh(toRemove) = [];
end
%% filter by UID input
if ~isempty(UID)
        [toRemove] = ~ismember(spikes.UID,UID);
    spikes.UID(toRemove) = [];
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.times{r} = [];
        end
    end
    spikes.times = removeEmptyCells(spikes.times);
    spikes.cluID(toRemove) = [];
    spikes.shankID(toRemove) = [];
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.rawWaveform{r} = [];
        end
    end
    spikes.rawWaveform = removeEmptyCells(spikes.rawWaveform);
    spikes.maxWaveformCh(toRemove) = [];
end

%% save to buzcode format
if saveMat
    save([basepath filesep meta.FileName '.spikes.cellinfo.mat'],'spikes')
end






