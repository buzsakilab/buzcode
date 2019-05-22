function nwb = addUnitsInfo(xml,nwb)
    % Add the units info from the standard Buzcode format: *.spikes.cellinfo.mat
    % This code was tested on the Buzcode tutorial dataset: 20170505_396um_0um_merge
    % Konstantinos Nasiotis 2019
    
    %%
    cellinfoFiles = dir([xml.folder_path filesep '*spikes.cellinfo.mat']);
    
    if length(cellinfoFiles)>1
        disp('More than one spikes.cellinfo file is present. Exiting without adding spiking information')
        return
    elseif isempty(cellinfoFiles)
        disp('No spikes.cellinfo.mat files present in the selected directory.')
        return
    end
    
    % If there a single spikes.cellinfo.mat file, load it
    load(cellinfoFiles.name)

    %% Add spikes info
    % Serialize spiketimes and cluIDs
    spike_times       = [];
    spike_times_index = [];

    current_index = 0;
    for iNeuron = 1:length(spikes.UID)
        spike_times = [spike_times ; spikes.times{iNeuron}];
        spike_times_index = [spike_times_index; int64(length(spikes.times{iNeuron})+current_index)];
        current_index = spike_times_index(end);
    end

    electrode_group = [];
    shank_that_neurons_belongs_to = zeros(length(spikes.UID),1);
    for iNeuron = 1:length(spikes.UID)
        electrode_group = [electrode_group; types.untyped.ObjectView(['/general/extracellular_ephys/' ['shank' num2str(spikes.shankID(iNeuron))]])];
    end

    electrode_group = types.core.VectorData('data', electrode_group, 'description','the electrode group that each spike unit came from');

    % Initialize the fields needed
    spike_times       = types.core.VectorData        ('data', spike_times, 'description', 'the spike times for each unit');
    spike_times_index = types.core.VectorIndex       ('data', spike_times_index, 'target', types.untyped.ObjectView('/units/spike_times')); % The ObjectView links the indices to the spike times
    id                = types.core.ElementIdentifiers('data', spikes.UID');

    
    waveform_mean = zeros(length(spikes.UID),length(spikes.rawWaveform{1}));
    for iNeuron = 1:length(spikes.UID)
        waveform_mean(iNeuron,:) = spikes.rawWaveform{iNeuron};
    end
    waveform_mean = types.core.VectorData('data', waveform_mean, 'description', 'The mean Waveform for each unit. nNeurons x nSamples');
    
    

    %% Fill the units fields
    nwb.units = types.core.Units( ...
        'electrode_group', electrode_group, 'electrodes', [], 'electrodes_index', [], 'obs_intervals', [], 'obs_intervals_index', [], ...
        'spike_times', spike_times, 'spike_times_index', spike_times_index, 'waveform_mean', waveform_mean, 'waveform_sd', [], ...
        'colnames', {'shank_id'; 'spike_times'; 'electrode_group'; 'cell_type'; 'global_id'; 'max_electrode'}, ...
        'description', 'Generated from Buzcode2NWB - addUnitsInfo', 'id', id, 'vectorindex', []);         

    %% Extra Unit Info
    nwb.units.vectordata.set('shankID'      , types.core.VectorData('description', 'Which shank each unit belongs to', 'data', spikes.shankID));
    nwb.units.vectordata.set('cluID'        , types.core.VectorData('description', 'cluster ID', 'data', spikes.cluID));
    nwb.units.vectordata.set('maxWaveformCh', types.core.VectorData('description', 'The electrode where each unit showed maximum Waveform', 'data', spikes.maxWaveformCh));
    nwb.units.vectordata.set('region'       , types.core.VectorData('description', 'The region where each neuron belongs to', 'data', spikes.region));

    
    disp('Spikes info added..')
    
    
end