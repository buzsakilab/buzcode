function spikes = bz_GetSpikes(varargin)
% bz_getSpikes - Get spike timestamps.
%       if loading from clu/res/fet/spk files - must be formatted as:
%       baseName.clu.shankNum
%
% USAGE
%
%    spikes = bz_getSpikes(varargin)
% 
% INPUTS
%
%    spikeGroups     -vector subset of shank IDs to load (Default: all)
%    region          -string region ID to load neurons from specific region
%                     (requires sessionInfo file or units->structures in xml)
%    UID             -vector subset of UID's to load 
%    basepath        -path to recording (where .dat/.clu/etc files are)
%    getWaveforms    -logical (default=true) to load mean of raw waveform data
%    forceReload     -logical (default=false) to force loading from
%                     res/clu/spk files
%    onlyLoad        -[shankID cluID] pairs to EXCLUSIVELY LOAD from 
%                       clu/res/fet to spikes.cellinfo.mat file
%    saveMat         -logical (default=false) to save in buzcode format
%    noPrompts       -logical (default=false) to supress any user prompts
%    
% OUTPUTS
%
%    spikes - cellinfo struct with the following fields
%          .sessionName    -name of recording file
%          .UID            -unique identifier for each neuron in a recording
%          .times          -cell array of timestamps (seconds) for each neuron
%          .spindices      -sorted vector of [spiketime UID], useful for 
%                           input to some functions and plotting rasters
%          .region         -region ID for each neuron (especially important large scale, high density probes)
%          .shankID        -shank ID that each neuron was recorded on
%          .maxWaveformCh  -channel # with largest amplitude spike for each neuron
%          .rawWaveform    -average waveform on maxWaveformCh (from raw .dat)
%          .cluID          -cluster ID, NOT UNIQUE ACROSS SHANKS
%          .numcells       -number of cells/UIDs
%           
% NOTES
%
% This function can be used in several ways to load spiking data.
% Specifically, it loads spiketimes for individual neurons and other
% sessionInfodata that describes each neuron.  Spiketimes can be loaded using the
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
%% Deal With Inputs 
spikeGroupsValidation = @(x) assert(isnumeric(x) || strcmp(x,'all'),...
    'spikeGroups must be numeric or "all"');

p = inputParser;
addParameter(p,'spikeGroups','all',spikeGroupsValidation);
addParameter(p,'region','',@isstr); % won't work without sessionInfodata 
addParameter(p,'UID',[],@isvector);
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'getWaveforms',true)
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'onlyLoad',[]);

parse(p,varargin{:})

spikeGroups = p.Results.spikeGroups;
region = p.Results.region;
UID = p.Results.UID;
basepath = p.Results.basepath;
getWaveforms = p.Results.getWaveforms;
forceReload = p.Results.forceReload;
saveMat = p.Results.saveMat;
noPrompts = p.Results.noPrompts;
onlyLoad = p.Results.onlyLoad;


[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', noPrompts);
baseName = bz_BasenameFromBasepath(basepath);


spikes.samplingRate = sessionInfo.rates.wideband;
nChannels = sessionInfo.nChannels;

cellinfofile = [basepath filesep sessionInfo.FileName '.spikes.cellinfo.mat'];
%% if the cellinfo file exist and we don't want to re-load files
if exist(cellinfofile,'file') && forceReload == false
    disp('loading spikes from cellinfo file..')
    load(cellinfofile)
    %Check that the spikes structure fits cellinfo requirements
    [iscellinfo] = bz_isCellInfo(spikes);
    switch iscellinfo
        case false
            warning(['The spikes structure in baseName.spikes.cellinfo.mat ',...
                'does not fit buzcode standards. Sad.'])
    end
    
    %If regions have been added since creation... add them
    if ~isfield(spikes,'region') & isfield(sessionInfo,'region')
        if ~isfield(spikes,'numcells')
            spikes.numcells = length(spikes.UID);
        end
        for cc = 1:spikes.numcells
            spikes.region{cc} = sessionInfo.region{spikes.maxWaveformCh(cc)==sessionInfo.channels};
        end
        
        if saveMat
            save(cellinfofile,'spikes')
        end
    end
    
else % do the below then filter by inputs... (Load from clu/res/fet)
    
    if ~noPrompts & saveMat == 0 %Inform the user that they should save a file for later
        savebutton = questdlg(['Would you like to save your spikes in ',...
            sessionInfo.FileName,'.spikes.cellinfo.mat?  ',...
            'This will save significant load time later.']);
        if strcmp(savebutton,'Yes'); saveMat = true; end
    end
    
disp('loading spikes from clu/res/spk files..')
% find res/clu/fet/spk files here
cluFiles = dir([basepath filesep '*.clu*']);  
resFiles = dir([basepath filesep '*.res.*']);
if any(getWaveforms)
    spkFiles = dir([basepath filesep '*.spk*']);
end


% remove *temp*, *autosave*, and *.clu.str files/directories
tempFiles = zeros(length(cluFiles),1);
for i = 1:length(cluFiles) 
    dummy = strsplit(cluFiles(i).name, '.'); % Check whether the component after the last dot is a number or not. If not, exclude the file/dir. 
    if ~isempty(findstr('temp',cluFiles(i).name)) | ~isempty(findstr('autosave',cluFiles(i).name)) | ...
            isempty(str2num(dummy{length(dummy)})) | find(contains(dummy, 'clu')) ~= length(dummy)-1  | ...
            ~strcmp(dummy{1},baseName)
        tempFiles(i) = 1;
    end
end
cluFiles(tempFiles==1)=[];
tempFiles = zeros(length(resFiles),1);
for i = 1:length(resFiles)
    dummy = strsplit(resFiles(i).name, '.');
    if ~isempty(findstr('temp',resFiles(i).name)) | ~isempty(findstr('autosave',resFiles(i).name)) | ...
            ~strcmp(dummy{1},baseName)
        tempFiles(i) = 1;
    end
end
if any(getWaveforms)
    resFiles(tempFiles==1)=[];
    tempFiles = zeros(length(spkFiles),1);
    for i = 1:length(spkFiles)
        dummy = strsplit(spkFiles(i).name, '.');
        if ~isempty(findstr('temp',spkFiles(i).name)) | ~isempty(findstr('autosave',spkFiles(i).name)) | ...
            ~strcmp(dummy{1},baseName)
            tempFiles(i) = 1;
        end
    end
    spkFiles(tempFiles==1)=[];
end

if isempty(cluFiles)
    disp('no clu files found...')
    spikes = [];
    return
end


% ensures we load in sequential order (forces compatibility with FMAT
% ordering)
for i = 1:length(cluFiles)
    temp = strsplit(cluFiles(i).name,'.');
    shanks(i) = str2num(temp{length(temp)});
end
[shanks ind] = sort(shanks);
cluFiles = cluFiles(ind); %Bug here if there are any files x.clu.x that are not your desired clus
resFiles = resFiles(ind);
if any(getWaveforms)
    spkFiles = spkFiles(ind);
end

% check if there are matching #'s of files
if length(cluFiles) ~= length(resFiles) & length(cluFiles) ~= length(spkFiles)
    error('found an incorrect number of res/clu/spk files...')
end

% use the .clu files to get spike ID's and generate UID and spikeGroup
% use the .res files to get spike times
count = 1;

if isempty(sessionInfo.spikeGroups.groups)
    sessionInfo.spikeGroups = sessionInfo.AnatGrps;
end
for i=1:length(cluFiles) 
    disp(['working on ' cluFiles(i).name])
    
    temp = strsplit(cluFiles(i).name,'.');
    shankID = str2num(temp{length(temp)}); %shankID is the spikegroup number
    clu = load(fullfile(basepath,cluFiles(i).name));
    clu = clu(2:end); % toss the first sample to match res/spk files
    res = load(fullfile(basepath,resFiles(i).name));
    spkGrpChans = sessionInfo.spikeGroups.groups{shankID}; % we'll eventually want to replace these two lines
    
    if any(getWaveforms) && sum(clu)>0 %bug fix if no clusters 
        nSamples = sessionInfo.spikeGroups.nSamples(shankID);
        % load waveforms
        chansPerSpikeGrp = length(sessionInfo.spikeGroups.groups{shankID});
        fid = fopen(fullfile(basepath,spkFiles(i).name),'r');
        wav = fread(fid,[1 inf],'int16=>int16');
        try %bug in some spk files... wrong number of samples?
            wav = reshape(wav,chansPerSpikeGrp,nSamples,[]);
        catch
            if strcmp(getWaveforms,'force')
                wav = nan(chansPerSpikeGrp,nSamples,length(clu));
                display([spkFiles(i).name,' error.'])
            else
            error(['something is wrong with ',spkFiles(i).name,...
                ' Use ''getWaveforms'', false to skip waveforms or ',...
                '''getWaveforms'', ''force'' to write nans on bad shanks.'])
            end
        end
        wav = permute(wav,[3 1 2]);
    end
    
    cells  = unique(clu);
    % remove MUA and NOISE clusters...
    cells(cells==0) = [];
    cells(cells==1) = [];  % consider adding MUA as another input argument...?
    
    for c = 1:length(cells)
       spikes.UID(count) = count; % this only works if all shanks are loaded... how do we optimize this?
       ind = find(clu == cells(c));
       spikes.times{count} = res(ind) ./ spikes.samplingRate;
       spikes.shankID(count) = shankID;
       spikes.cluID(count) = cells(c);

       %Waveforms    
       if any(getWaveforms)
           wvforms = squeeze(mean(wav(ind,:,:)))-mean(mean(mean(wav(ind,:,:)))); % mean subtract to account for slower (theta) trends
           if prod(size(wvforms))==length(wvforms)%in single-channel groups wvforms will squeeze too much and will have amplitude on D1 rather than D2
               wvforms = wvforms';%fix here
           end
           for t = 1:size(wvforms,1)
              [a(t) b(t)] = max(abs(wvforms(t,:))); 
           end
           [aa bb] = max(a,[],2);
           spikes.rawWaveform{count} = wvforms(bb,:);
           spikes.maxWaveformCh(count) = spkGrpChans(bb);  
           %Regions (needs waveform peak)
           if isfield(sessionInfo,'region') %if there is regions field in your metadata
                spikes.region{count} = sessionInfo.region{find(spkGrpChans(bb)==sessionInfo.channels)};
           elseif isfield(sessionInfo,'Units') %if no regions, but unit region from xml via Loadparamteres
                %Find the xml Unit that matches group/cluster
                unitnum = cellfun(@(X,Y) X==spikes.shankID(count) && Y==spikes.cluID(count),...
                    {sessionInfo.Units(:).spikegroup},{sessionInfo.Units(:).cluster});
                if sum(unitnum) == 0
                    display(['xml Missing Unit - spikegroup: ',...
                        num2str(spikes.shankID(count)),' cluster: ',...
                        num2str(spikes.cluID(count))])
                    spikes.region{count} = 'missingxml';
                else %possible future bug: two xml units with same group/clu...              
                    spikes.region{count} = sessionInfo.Units(unitnum).structure;
                end
           end
           clear a aa b bb
       end
       
       count = count + 1;
    end
end


if ~isempty(onlyLoad)
    toRemove = true(size(spikes.UID));
    for cc = 1:size(onlyLoad,1)
        whichUID = ismember(spikes.shankID,onlyLoad(cc,1)) & ismember(spikes.cluID,onlyLoad(cc,2));
        toRemove(whichUID) = false;
        if ~any(whichUID)
            display(['No unit with shankID:',num2str(onlyLoad(cc,1)),...
                ' cluID:',num2str(onlyLoad(cc,2))])
        end
    end
    spikes = removeCells(toRemove,spikes,getWaveforms);
end


spikes.sessionName = sessionInfo.FileName;
end

%% save to buzcode format (before exclusions)
if saveMat
    save(cellinfofile,'spikes')
end


%% EXCLUSIONS %%

%filter by spikeGroups input
if ~strcmp(spikeGroups,'all')
    [toRemove] = ~ismember(spikes.shankID,spikeGroups);
    spikes = removeCells(toRemove,spikes,getWaveforms);
end

%filter by region input
if ~isempty(region)
    if ~isfield(spikes,'region') %if no region information in metadata
        error(['You selected to load cells from region "',region,...
            '", but there is no region information in your sessionInfo'])
    end
    
  toRemove = ~ismember(spikes.region,region);
    if sum(toRemove)==length(spikes.UID) %if no cells from selected region
        warning(['You selected to load cells from region "',region,...
            '", but none of your cells are from that region'])
    end
    
    spikes = removeCells(toRemove,spikes,getWaveforms);
end

%filter by UID input
if ~isempty(UID)
	[toRemove] = ~ismember(spikes.UID,UID);
    spikes = removeCells(toRemove,spikes,getWaveforms);   
end

%% Generate spindices matrics
spikes.numcells = length(spikes.UID);
for cc = 1:spikes.numcells
    groups{cc}=spikes.UID(cc).*ones(size(spikes.times{cc}));
end
if spikes.numcells>0
    alltimes = cat(1,spikes.times{:}); groups = cat(1,groups{:}); %from cell to array
    [alltimes,sortidx] = sort(alltimes); groups = groups(sortidx); %sort both
    spikes.spindices = [alltimes groups];
end

%% Check if any cells made it through selection
if isempty(spikes.times) | spikes.numcells == 0
    spikes = [];
end



end


%%
function spikes = removeCells(toRemove,spikes,getWaveforms)
%Function to remove cells from the structure. toRemove is the INDEX of
%the UID in spikes.UID

    spikes.UID(toRemove) = [];
    spikes.times(toRemove) = [];
    spikes.region(toRemove) = [];
    spikes.cluID(toRemove) = [];
    spikes.shankID(toRemove) = [];
    
    if any(getWaveforms)
        spikes.rawWaveform(toRemove) = [];
        spikes.maxWaveformCh(toRemove) = [];
    end
    
end





