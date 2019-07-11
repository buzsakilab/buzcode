 function nwb = addUnitsInfo_Neuroscope(xml, nwb)
    %% Add the units info (copied from bz_GetSpikes)
    % Adds unit info in: nwb.units

    % This code takes the unit information from the .clu, .res, .spk files
    
    
    %% Get unit information
    getWaveforms = 1; % Set this to true if you want to add waveforms on the NWB file


    spikes.samplingRate = str2double(xml.acquisitionSystem.samplingRate.Text);


    disp('loading spikes from clu/res/spk files..')
    % find res/clu/fet/spk files here
    cluFiles = dir([xml.folder_path filesep '*.clu*']);  
    resFiles = dir([xml.folder_path filesep '*.res*']);
    if any(getWaveforms)
        spkFiles = dir([xml.folder_path filesep '*.spk*']);
    end

    % remove *temp*, *autosave*, and *.clu.str files/directories
    tempFiles = zeros(length(cluFiles),1);
    for i = 1:length(cluFiles) 
        dummy = strsplit(cluFiles(i).name, '.'); % Check whether the component after the last dot is a number or not. If not, exclude the file/dir. 
        if ~isempty(findstr('temp',cluFiles(i).name)) | ~isempty(findstr('autosave',cluFiles(i).name)) | isempty(str2num(dummy{length(dummy)})) | find(contains(dummy, 'clu')) ~= length(dummy)-1  
            tempFiles(i) = 1;
        end
    end
    cluFiles(tempFiles==1)=[];
    tempFiles = zeros(length(resFiles),1);
    for i = 1:length(resFiles)
        if ~isempty(findstr('temp',resFiles(i).name)) | ~isempty(findstr('autosave',resFiles(i).name))
            tempFiles(i) = 1;
        end
    end
    if any(getWaveforms)
        resFiles(tempFiles==1)=[];
        tempFiles = zeros(length(spkFiles),1);
        for i = 1:length(spkFiles)
            if ~isempty(findstr('temp',spkFiles(i).name)) | ~isempty(findstr('autosave',spkFiles(i).name))
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
    if length(cluFiles) ~= length(resFiles) && length(cluFiles) ~= length(spkFiles)
        error('found an incorrect number of res/clu/spk files...')
    end

    % use the .clu files to get spike ID's and generate UID and spikeGroup
    % use the .res files to get spike times
    count = 1;

    ecephys = types.core.ProcessingModule;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is copied from the ElectrodesInfo
    nShanks = length(xml.spikeDetection.channelGroups.group);
    groups = xml.spikeDetection.channelGroups.group; % Use this for simplicity
    all_shank_channels = cell(nShanks,1); % This will hold the channel numbers that belong in each shank
    shank = [];
    group_object_view = [];

    for iGroup = 1:nShanks
    % Get all_shank_channls again  for iGroup = 1:nShanks
        for iChannel = 1:length(groups{iGroup}.channels.channel)
            all_shank_channels{iGroup} = [all_shank_channels{iGroup} str2double(groups{iGroup}.channels.channel{iChannel}.Text)];
            shank = [shank iGroup];
            group_object_view = [group_object_view; types.untyped.ObjectView(['/general/extracellular_ephys/' ['shank' num2str(iGroup)]])];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    for iShank=1:length(cluFiles) 
        disp(['working on ' cluFiles(iShank).name])

        temp = strsplit(cluFiles(iShank).name,'.');
        shankID = str2num(temp{length(temp)}); %shankID is the spikegroup number
        clu = load(fullfile(xml.folder_path,cluFiles(iShank).name));
        clu = clu(2:end); % toss the first sample to match res/spk files
        res = load(fullfile(xml.folder_path,resFiles(iShank).name));
        spkGrpChans = all_shank_channels{iShank};

        if any(getWaveforms) && sum(clu)>0 %bug fix if no clusters 
            nSamples = str2double(xml.spikeDetection.channelGroups.group{iShank}.nSamples.Text);
            % load waveforms
            chansPerSpikeGrp = length(all_shank_channels{iShank});
            fid = fopen(fullfile(xml.folder_path,spkFiles(iShank).name),'r');
            wav = fread(fid,[1 inf],'int16=>int16');
            try %bug in some spk files... wrong number of samples?
                wav = reshape(wav,chansPerSpikeGrp,nSamples,[]);
            catch
                if strcmp(getWaveforms,'force')
                    wav = nan(chansPerSpikeGrp,nSamples,length(clu));
                    display([spkFiles(iShank).name,' error.'])
                else
                error(['something is wrong with ',spkFiles(iShank).name,...
                    ' Use ''getWaveforms'', false to skip waveforms or ',...
                    '''getWaveforms'', ''force'' to write nans on bad shanks.'])
                end
            end
            wav = permute(wav,[3 1 2]);

            %% Get the DynamicTableRegion field for each shank

            % First check if the electrodes field has been filled
            if isempty(nwb.general_extracellular_ephys_electrodes)
                nwb = Neuroscope2NWB.getElectrodeInfo(xml, nwb);
            end

            electrodes_field = types.core.DynamicTableRegion('table',types.untyped.ObjectView('/general/extracellular_ephys/electrodes'),'description',['shank' num2str(iShank) ' region'],'data',nwb.general_extracellular_ephys_electrodes.id.data(find(shank == iShank)'));
            SpikeEventSeries = types.core.SpikeEventSeries('data', wav, 'electrodes', electrodes_field, 'timestamps', res./ spikes.samplingRate);

            %% This section assigns the spike-waveforms in the .NWB
            ecephys.nwbdatainterface.set(['SpikeEventSeries' num2str(iShank)],SpikeEventSeries);

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
               spikes.maxWaveformCh(count) = spkGrpChans(bb);  % Use this in Brainstorm
    %            %Regions (needs waveform peak)
    %            if isfield(xml,'region') %if there is regions field in your metadata
    %                 spikes.region{count} = 'unknown';
    %            elseif isfield(xml,'units') %if no regions, but unit region from xml via Loadparamteres
    %                 %Find the xml Unit that matches group/cluster
    %                 unitnum = cellfun(@(X,Y) X==spikes.shankID(count) && Y==spikes.cluID(count),...
    %                     {sessionInfo.Units(:).spikegroup},{sessionInfo.Units(:).cluster});
    %                 if sum(unitnum) == 0
    %                     display(['xml Missing Unit - spikegroup: ',...
    %                         num2str(spikes.shankID(count)),' cluster: ',...
    %                         num2str(spikes.cluID(count))])
    %                     spikes.region{count} = 'missingxml';
    %                 else %possible future bug: two xml units with same group/clu...              
    %                     spikes.region{count} = sessionInfo.Units(unitnum).structure;
    %                 end
    %            end
               clear a aa b bb
           end

           count = count + 1;

        end

        ecephys.description = 'intermediate data from extracellular electrophysiology recordings, e.g., LFP';
        nwb.processing.set('ecephys', ecephys);
    end


    % Serialize spiketimes and cluIDs
    spike_times       = [];
    spike_times_index = [];

    current_index = 0;
    for iNeuron = 1:length(spikes.UID)
        spike_times = [spike_times ; spikes.times{iNeuron}];
        spike_times_index = [spike_times_index; int64(length(spikes.times{iNeuron})+current_index)];
        current_index = spike_times_index(end);
    end


    % electrode_group - Assigns the group_object_view that was defined above at
    % the electrodes, to specific neurons - I need to find how each neuron is
    % assigned to a shank
    electrode_group = [];
    shank_that_neurons_belongs_to = zeros(length(spikes.UID),1);
    for iNeuron = 1:length(spikes.UID)
        shank_that_neurons_belongs_to(iNeuron) = str2double(xml.units.unit{iNeuron}.group.Text);
        first_electrode_in_shank = find(shank == shank_that_neurons_belongs_to(iNeuron));
        first_electrode_in_shank = first_electrode_in_shank(1);
        electrode_group = [electrode_group; group_object_view(first_electrode_in_shank)];
    end

    electrode_group = types.core.VectorData('data', electrode_group, 'description','the electrode group that each spike unit came from');

    % Initialize the fields needed
    spike_times       = types.core.VectorData        ('data', spike_times, 'description', 'the spike times for each unit');
    spike_times_index = types.core.VectorIndex       ('data', spike_times_index, 'target', types.untyped.ObjectView('/units/spike_times')); % The ObjectView links the indices to the spike times
    id                = types.core.ElementIdentifiers('data', [0:length(xml.units.unit)-1]');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THIS GAVE AN ERROR WHEN ASSIGNING CELL ARRAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    waveform_mean = types.core.VectorData('data', spikes.rawWaveform, 'description', 'The mean Waveform for each unit');
    waveform_mean = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% Fill the units fields
    nwb.units = types.core.Units( ...
        'electrode_group', electrode_group, 'electrodes', [], 'electrodes_index', [], 'obs_intervals', [], 'obs_intervals_index', [], ...
        'spike_times', spike_times, 'spike_times_index', spike_times_index, 'waveform_mean', waveform_mean, 'waveform_sd', [], ...
        'colnames', {'shank_id'; 'spike_times'; 'electrode_group'; 'cell_type'; 'global_id'; 'max_electrode'}, ...
        'description', 'Generated from Neuroscope2NWB', 'id', id, 'vectorindex', []);         

    %% Extra Unit Info
    % FOR THE VECTORDATA, IDEALLY I NEED FILE: DG_all_6__UnitFeatureSummary_add (ACCORDING TO BEN'S CONVERTER - THIS HOLDS INFO ABOUT THE CELL_TYPE, GLOBAL_ID)
    nwb.units.vectordata.set('cluID',         types.core.VectorData('description', 'cluster ID', 'data', spikes.cluID));
    nwb.units.vectordata.set('maxWaveformCh', types.core.VectorData('description', 'The electrode where each unit showed maximum Waveform', 'data', spikes.maxWaveformCh));

    disp('Spikes info added..')
end