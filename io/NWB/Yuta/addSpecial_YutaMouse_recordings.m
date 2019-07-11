function nwb = addSpecial_YutaMouse_recordings(xml,nwb)
    %% Add raw recordings in: nwb.acquisition

    % values taken from Yuta's spreadsheet
    % HOW ABOUT POSITION0 - POSITION1 CHANNELS???

    hdr.nChannels = str2double(xml.acquisitionSystem.nChannels.Text);
    hdr.sRateLfp  = str2double(xml.fieldPotentials.lfpSamplingRate.Text);

    lfpFile = dir([xml.folder_path filesep '*.eeg']);

    special_electrode_labels  = {'ch_wait','ch_arm','ch_solL','ch_solR','ch_dig1','ch_dig2','ch_entL','ch_entR','ch_SsolL','ch_SsolR'};
    special_electrode_indices = [79,78,76,77,65,68,72,71,73,70]; 

    for iSpecialElectrode = 1:length(special_electrode_labels)
        special_Electrode_data = bz_LoadBinary(fullfile(lfpFile.folder,lfpFile.name), 'duration',Inf, 'frequency',hdr.sRateLfp,'nchannels',hdr.nChannels, 'channels', special_electrode_indices(iSpecialElectrode));
        single_Electrode = types.core.TimeSeries('description','environmental electrode recorded inline with neural data','data',special_Electrode_data,'starting_time', 0, 'starting_time_rate', hdr.sRateLfp, 'data_unit','V');
        nwb.acquisition.set(special_electrode_labels{iSpecialElectrode}, single_Electrode);
    end

    disp('Special YutaMouse channel info added..')
end