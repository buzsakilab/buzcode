function nwb = addElectrophysiology(xml,nwb)
    %% Adds electrophysiology signals in: nwb.processing.get('ecephys').nwbdatainterface.get('LFP')
    % This code was tested on the Buzcode tutorial dataset: 20170505_396um_0um_merge
    % and the YutaMouse41-150903
    % Konstantinos Nasiotis 2019

    %%
    [ff basename] = fileparts(xml.folder_path);


    % bz_LoadBinary
    lfpFile = dir([xml.folder_path filesep basename '*.lfp']);

    if length(lfpFile)>1
        disp('More than one .eeg files are present here. No Electrophysiology signals were added')
        return
    elseif length(lfpFile)==0
        lfpFile = dir([xml.folder_path filesep basename '*.eeg']);
        if length(lfpFile)>1
           disp('More than one .lfp files are present here. No Electrophysiology signals were added')
           return
        elseif length(lfpFile)==0
           disp('No .eeg or .lfp files are present in the selected directory. No Electrophysiology signals were added')
           return
        end
    end


    % Get the samples number, based on the size of the file
    % Check for the precision that samples are saved first

    hdr.nBits     = str2double(xml.acquisitionSystem.nBits.Text);
    hdr.nChannels = str2double(xml.acquisitionSystem.nChannels.Text);
    hdr.sRateOrig = str2double(xml.acquisitionSystem.samplingRate.Text);
    hdr.Gain      = str2double(xml.acquisitionSystem.amplification.Text);
    hdr.sRateLfp  = str2double(xml.fieldPotentials.lfpSamplingRate.Text);


    % Get data type
    switch lower(hdr.nBits)
        case 16;
            hdr.byteSize   = 2;
            hdr.byteFormat = 'int16';
        case 32;
            hdr.byteSize   = 4;
            hdr.byteFormat = 'int32';
    end
    % Guess the number of time points based on the file size
    dirInfo = dir(fullfile(lfpFile.folder,lfpFile.name));
    hdr.nSamples = floor(dirInfo.bytes ./ (hdr.nChannels * hdr.byteSize));

    lfp_data = bz_LoadBinary(fullfile(lfpFile.folder,lfpFile.name), 'duration',Inf, 'frequency',hdr.sRateLfp,'nchannels',hdr.nChannels, 'channels', [1:hdr.nChannels]); % nSamples x 64

    % If the electrode Information has not already been filled, 
    % do it now
    if isempty(nwb.general_extracellular_ephys_electrodes)
        nwb = Neuroscope2NWB.getElectrodeInfo(xml,nwb);
    end


    electrodes_field = types.core.DynamicTableRegion('table',types.untyped.ObjectView('/general/extracellular_ephys/electrodes'),'description','electrode table reference','data',nwb.general_extracellular_ephys_electrodes.id.data);

    lfp = types.core.ElectricalSeries('data', lfp_data', 'electrodes',electrodes_field, 'description', 'lfp signal for all shank electrodes', 'starting_time', 0, 'starting_time_rate', hdr.sRateLfp);
    % I TRANSPOSED THE MATRIX HERE TO BE COMPATIBLE WITH THE FUNCTION BZ_GET_LFP

    LFP = types.core.LFP;
    LFP.electricalseries.set('lfp',lfp); 

    % Check if the ecephys field is already created           
    if isempty(keys(nwb.processing))
        ecephys = types.core.ProcessingModule('description', '');
        ecephys.description = 'intermediate data from extracellular electrophysiology recordings, e.g., LFP';
        nwb.processing.set('ecephys', ecephys);
    else
        if ~ismember(keys(nwb.processing),'ecephys')
            ecephys = types.core.ProcessingModule('description', '');
            ecephys.description = 'intermediate data from extracellular electrophysiology recordings, e.g., LFP';
            nwb.processing.set('ecephys', ecephys);
        end
    end

    nwb.processing.get('ecephys').nwbdatainterface.set('LFP', LFP);
    disp('Electrophysiological signals added..')
end