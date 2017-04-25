function spikes = bz_getSpikes(varargin)
% bz_getSpikes - Get spike timestamps.
%
% USAGE
%
%    spikes = bz_getSpikes(units,<options>)
% 
% INPUTS
%
%    spikeGroups     -vector subset of shank IDs to load
%    region          -string region ID to load neurons from specific region (requires metadata file)
%    UID             -vector subset of UID's to load 
%    basepath        -path to recording (where .dat/.clu/etc files are)
%    getwav    -logical (default=true) to load waveform data
%    forceReload     -logical (default=false) to force loading from
%                     res/clu/spk files
%    
% OUTPUTS
%
%    spikes - cellinfo struct with the following fields
%          .sessionName    -name of recording file
%          .UID            -unique identifier for each neuron in a recording
%          .times          -cell array of timestamps (seconds) for each neuron
%          .region         -region ID for each neuron (especially important large scale, high density probes)
%          .spikeGroup     -shank ID that each neuron was recorded on
%          .maxWaveformCh  -channel # with largest amplitude spike for each neuron
%          .meanWaveform   -average waveform on maxWaveformCh
%          .cluID          -cluster ID, NOT UNIQUE ACROSS SHANKS
%           
% NOTES
%
% written by David Tingley, 2017

p = inputParser;
addParameter(p,'spikeGroups',[],@isvector);
addParameter(p,'region','',@isstr); % won't work without metadata 
addParameter(p,'UID',[],@isvector);
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'getwav',true,@islogical)
addParameter(p,'forceReload',false,@islogical);


parse(p,varargin{:})

spikeGroups = p.Results.spikeGroups;
region = p.Results.region;
UID = p.Results.UID;
basepath = p.Results.basepath;
forceReload = p.Results.forceReload;

% get meta info about recording, for now we'll use the xml.
meta = LoadParameters(basepath);
samplingRate = meta.rates.wideband;
nChannels = meta.nChannels;


% to use session metadata, ADD IT HERE AS 'meta' and not below.

% samplingRate = meta.rates.wdieband;
% nChannels = meta.nChannels;


%% if the cellinfo file exist and we don't want to re-load files
if exist([basepath filesep 'spikes.cellinfo.mat']) & forceReload == false
    load([basepath filesep 'spikes.cellinfo.mat'])
    return
end

% find res/clu/fet/spk files here
cluFiles = dir([basepath filesep '*.clu*']);
resFiles = dir([basepath filesep '*.res*']);
spkFiles = dir([basepath filesep '*.spk*']);

% check if there are matching #'s of files
if length(cluFiles) ~= length(resFiles) & length(cluFiles) ~= length(spkFiles)
    error('found an incorrect number of res/clu/spk files...')
end

% use the .clu files to get spike ID's and generate UID and spikeGroup
% use the .res files to get spike times
count = 1;

for i=1:length(cluFiles)
    clu = load(cluFiles(i).name);
    clu = clu(2:end); % toss the first sample to match res/spk files
    res = load(resFiles(i).name);
    nSamples = meta.spikeGroups.nSamples(i);
    
    cells  = unique(clu);
    for c = 1:length(cells)
       spikes.UID(count) = count; % this only works if all shanks are loaded... how do we optimize this?
       
       ind = find(clu == cells(c));
       spikes.times{count} = res(ind) ./ samplingRate;
       spikes.spikeGroup(count) = i;
       spikes.cluID(count) = cells(c);
       
       % load waveforms
       
       wav = LoadSpikeWaveforms(spkFiles(i).name,nChannels,nSamples,ind(1:100));
    
       
       
       count = count + 1;
    end
end

spikes.sessionName = meta.FileName;


%% filter by spikeGroups input

%% filter by region input

%% filter by UID input












