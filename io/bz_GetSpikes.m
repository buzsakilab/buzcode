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
%    verbose         -logical (default=false)
%    keepCluWave     -logical (default=false) to keep waveform from .spk files
%                       as in previous bz_getSpikes functions (before 2019)
%                       If false (default), it instead uses waveforms directly from
%                       .dat files.  
%    sortingMethod   - [], 'kilosort' or 'clu'. If [], tries to detect a
%                   kilosort folder or clu files. 
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
%          .filtWaveform   -average filtered waveform on maxWaveformCh
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
% added Phy loading by Manu Valero, 2019 (previos bz_LoadPhy)

% TO DO: Get waveforms by an independent function (ie getWaveform) that
% generates a waveform.cellinfo.mat file with all channels waves.
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
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'onlyLoad',[]);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'keepCluWave',false,@islogical);
addParameter(p,'sortingMethod',[],@isstr);

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
verbose = p.Results.verbose;
keepCluWave = p.Results.keepCluWave;
sortingMethod = p.Results.sortingMethod;

[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', noPrompts);
baseName = bz_BasenameFromBasepath(basepath);

spikes.samplingRate = sessionInfo.rates.wideband;
nChannels = sessionInfo.nChannels;

cellinfofile = [basepath filesep sessionInfo.FileName '.spikes.cellinfo.mat'];
datfile = [basepath filesep sessionInfo.FileName '.dat'];
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
        if isfield(spikes,'maxWaveformCh')
            for cc = 1:spikes.numcells
                spikes.region{cc} = sessionInfo.region{spikes.maxWaveformCh(cc)==sessionInfo.channels};
            end
        end
        
        if saveMat
            save(cellinfofile,'spikes')
        end
    end
    
    if ~noPrompts & saveMat == 0 %Inform the user that they should save a file for later
        savebutton = questdlg(['Would you like to save your spikes in ',...
            sessionInfo.FileName,'.spikes.cellinfo.mat?  ',...
            'This will save significant load time later.']);
        if strcmp(savebutton,'Yes'); 
            saveMat = true; 
        end
    end
    
else
    % find res/clu/fet/spk files or kilosort folder here...
    cluFiles = dir([basepath filesep '*.clu*']);  
    resFiles = dir([basepath filesep '*.res.*']);
    if any(getWaveforms)
        spkFiles = dir([basepath filesep '*.spk*']);
    end
%     kilosort_path = dir([basepath filesep '*kilosort*']);
    %Finding kilosort folder, case insensitive - added by EFO
    temp = dir(basepath);
    aux  = strfind(lower({temp.name}),lower('kilosort'));
    fold_idx = find(cellfun(@(x) ~isempty(x),aux));
    %in case of multiple folders, it's going to take the first one in the
    %list
    kilosort_path = temp(fold_idx(1));
    
    if strcmpi(sortingMethod,'clu') || ~isempty(cluFiles) % LOADING FROM CLU/ RES
        fs = spikes.samplingRate; 
        disp('loading spikes from clu/res/spk files..')

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
        cellcount = 1;

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
               spikes.UID(cellcount) = cellcount; % this only works if all shanks are loaded... how do we optimize this?
               ind = find(clu == cells(c));
               spikes.times{cellcount} = res(ind) ./ spikes.samplingRate;
               spikes.shankID(cellcount) = shankID;
               spikes.cluID(cellcount) = cells(c);

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
                   spikes.rawWaveform{cellcount} = wvforms(bb,:);
                   spikes.maxWaveformCh(cellcount) = spkGrpChans(bb);  
                   clear a aa b bb
               end

               cellcount = cellcount + 1;
            end %end this cells in this clu            
        end%end this clu

        %%Regions: either using a) waveform peak, b) region of shank or c)%prior annotation in "structure" field of Units array
        if isfield(sessionInfo,'region') %if there is regions field in your metadata
            if isfield(spikes,'maxWaveformCh')% if max waveforms available, use region of peak channel
                for cellcount = 1:length(spikes.UID)%for each cell                    
                    spikes.region{cellcount} = sessionInfo.region{find(spikes.maxWaveformCh(cellcount)==sessionInfo.channels)};
                end
            else%if no max waveform, then use clu group and region of channels in that group to id region of unit
                for cellcount = 1:length(spikes.UID)%for each cell
                    tshank = spikes.shankID(cellcount);
                    if isfield(sessionInfo,'SpkGrps')
                        tshankchannels = sessionInfo.SpkGrps(tshank).Channels;
                    else
                        tshankchannels = sessionInfo.AnatGrps(tshank).Channels;
                    end
                    tshankchannelregions = sessionInfo.region(tshankchannels+1);
                    if prod(strcmp(tshankchannelregions{1},tshankchannelregions)) == 1 %if all the same
                        tregion = tshankchannelregions{1};
                    else
                        tregion = 'unknown'
                        disp(['For unit #' num2str(cellcount) ', cannot determine region.  The spike group it'' s found in has channels from different regions'])
                    end
                    spikes.region{cellcount} = tregion;
                end
            end
        elseif isfield(sessionInfo,'Units') %if no regions, but unit region from xml via Loadparamteres
            for cellcount = 1:length(spikes.UID)%for each cell
            %Find the xml Unit that matches group/cluster
                unitnum = cellfun(@(X,Y) X==spikes.shankID(cellcount) && Y==spikes.cluID(cellcount),...
                    {sessionInfo.Units(:).spikegroup},{sessionInfo.Units(:).cluster});
                if sum(unitnum) == 0
                    display(['xml Missing Unit - spikegroup: ',...
                        num2str(spikes.shankID(cellcount)),' cluster: ',...
                        num2str(spikes.cluID(cellcount))])
                    spikes.region{cellcount} = 'missingxml';
                else %possible future bug: two xml units with same group/clu...              
                    spikes.region{cellcount} = sessionInfo.Units(unitnum).structure;
                end
            end
        else
            disp('Unable to determine regions of units')
        end

        spikes.sessionName = sessionInfo.FileName;

    elseif strcmpi(sortingMethod, 'kilosort') || ~isempty(kilosort_path) % LOADING FROM KILOSORT

        disp('loading spikes from Kilosort/Phy format...')
        fs = spikes.samplingRate; 
        spike_cluster_index = readNPY(fullfile(kilosort_path.name, 'spike_clusters.npy'));
        spike_times = readNPY(fullfile(kilosort_path.name, 'spike_times.npy'));
        cluster_group = tdfread(fullfile(kilosort_path.name,'cluster_group.tsv'));
        
        %extracting shank information if it's kilosort2/phy2 output
        if exist(fullfile(kilosort_path.name,'cluster_info.tsv'),'file')
            cluster_info = tdfread(fullfile(kilosort_path.name,'cluster_info.tsv'));
            if isfield(cluster_info,'ch') %if it has the channel field
                clu_channels = cluster_info.ch;
                shanks = zeros(size(clu_channels));
                
                for s = 1:sessionInfo.spikeGroups.nGroups
                    temp1 = ismember(clu_channels,sessionInfo.spikeGroups.groups{s});
                    shanks(temp1) = s;
                end
            end
        else %otherwise try to load shanks.npy from kilosort1/phy1
            try
                shanks = readNPY(fullfile(kilosort_path.name, 'shanks.npy')); % done
            catch
                shanks = ones(size(cluster_group.cluster_id));
                warning('No shanks.npy file, assuming single shank!');
            end
        end
        spikes.sessionName = sessionInfo.FileName;
        jj = 1;
        for ii = 1:length(cluster_group.group)
            if strcmpi(strtrim(cluster_group.group(ii,:)),'good')
                ids = find(spike_cluster_index == cluster_group.cluster_id(ii)); % cluster id
                spikes.cluID(jj) = cluster_group.cluster_id(ii);
                spikes.UID(jj) = jj;
                spikes.times{jj} = double(spike_times(ids))/fs; % cluster time
                spikes.ts{jj} = double(spike_times(ids)); % cluster time
                cluster_id = find(cluster_group.cluster_id == spikes.cluID(jj));
                spikes.shankID(jj) = double(shanks(cluster_id));
                jj = jj + 1;
            end
        end

        if ~isfield(spikes,'region') && isfield(spikes,'maxWaveformCh') && isfield(sessionInfo,'region')
            for cc = 1:length(spikes.times)
                spikes.region{cc} = [sessionInfo.region{find(spikes.maxWaveformCh(cc)==sessionInfo.channels)} ''];
            end
        end
    else
        error('Unit format not recognized...');
    end

    % get waveforms from .dat file if possible, should be better than .spk
    % and as of 8/2020 is default to overwrite any .spk.  Will keep
    % .spk-based spikes if no .dat present however.
    if any(getWaveforms) && ~keepCluWave
        nPull = 1000;  % number of spikes to pull out
        wfWin = 0.008; % Larger size of waveform windows for filterning
        filtFreq = 500;
        hpFilt = designfilt('highpassiir','FilterOrder',3, 'PassbandFrequency',filtFreq,'PassbandRipple',0.1, 'SampleRate',fs);
        wfWin = round((wfWin * fs)/2);%in samples
        for ii = 1 : size(spikes.times,2)
            spkTmp = spikes.times{ii};
            if length(spkTmp) > nPull
                spkTmp = spkTmp(randperm(length(spkTmp)));
                spkTmp = spkTmp(1:nPull);
            end
            wf = [];
            for jj = 1 : length(spkTmp)
                if verbose
                    fprintf(' ** %3.i/%3.i for cluster %3.i/%3.i  \n',jj, length(spkTmp), ii, size(spikes.times,2));
                end
                %updated by EFO on 18/11/2020, bz_LoadBinary needs offset input in
                %samples and not in seconds
                wf = cat(3,wf,bz_LoadBinary([sessionInfo.session.name '.dat'],'offset',round(spkTmp(jj)*fs) - (wfWin),...
                    'samples',(wfWin * 2)+1,'frequency',sessionInfo.rates.wideband,'nChannels',sessionInfo.nChannels));
            end
            wf = mean(wf,3);
            if isfield(sessionInfo,'badchannels')
                wf(:,ismember(sessionInfo.channels,sessionInfo.badchannels))=0;
            end
            for jj = 1 : size(wf,2)          
                wfF(:,jj) = filtfilt(hpFilt,wf(:,jj) - mean(wf(:,jj)));
            end
            %updated by EFO on 18/11/2020 to only get the max channel on
            %the respective shank, avoiding getting waveform of a dead
            %channel
            shank_ch = sessionInfo.spikeGroups.groups{spikes.shankID(ii)}+1; %Channels are 0-based
            [~, maxCh] = max(abs(wfF(wfWin,shank_ch)));
            maxCh = shank_ch(maxCh);
            rawWaveform = detrend(wf(:,maxCh) - mean(wf(:,maxCh))); 
            filtWaveform = wfF(:,maxCh) - mean(wfF(:,maxCh));
            spikes.rawWaveform{ii} = rawWaveform(wfWin-(0.002*fs):wfWin+(0.002*fs)); % keep only +- 1ms of waveform
            spikes.filtWaveform{ii} = filtWaveform(wfWin-(0.002*fs):wfWin+(0.002*fs)); 
            spikes.maxWaveformCh(ii) = sessionInfo.channels(maxCh);
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

    %% save to buzcode format (before exclusions)
    if saveMat
        save(cellinfofile,'spikes')
    end

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
    spikes.shankID(toRemove) = [];
    if isfield(spikes,'cluID')
        spikes.cluID(toRemove) = [];
    elseif isfield(spikes,'UID_kilosort')
        spikes.UID_kilosort(toRemove) = [];
    end
    
    if any(getWaveforms)
        spikes.rawWaveform(toRemove) = [];
        spikes.maxWaveformCh(toRemove) = [];
        if isfield(spikes,'filtWaveform')
            spikes.filtWaveform(toRemove) = [];
        end
    end
end





