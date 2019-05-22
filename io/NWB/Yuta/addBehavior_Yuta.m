function nwb = addBehavior_Yuta(xml,nwb)

    %% Add behavioral data: nwb2.processing.get('behavior').nwbdatainterface
    % Check for YutaMouse behavioral files (*position*)

    behavioralFiles = dir(fullfile(xml.folder_path, '*position*'));

    if ~isempty(behavioralFiles)

        behavior_NWB = types.core.ProcessingModule;


        time_epochs = repmat(struct('label','','start_times',0,'stop_times',0),length(behavioralFiles),1);

        for iFile = 1:length(behavioralFiles)

            % The label of the behavior
            behavioral_Label = strsplit(behavioralFiles(iFile).name,'__');
            behavioral_Label = erase(behavioral_Label{2},'.mat');

            position_signals = load(fullfile(behavioralFiles(iFile).folder, behavioralFiles(iFile).name));

            % Some behavioral signals might have more than one signal in them
            field_names = fieldnames(position_signals);

            the_position_field_NWB = types.core.Position;

            for iField = 1:length(field_names)

                behavioral_timestamps = position_signals.(field_names{iField})(:,1);
                position_coordinates  = position_signals.(field_names{iField})(:,2:end);
                spatial_series = types.core.SpatialSeries('data', position_coordinates, 'timestamps', behavioral_timestamps, 'reference_frame', 'unknown', 'data_conversion', 1);

                the_position_field_NWB.spatialseries.set(field_names{iField}, spatial_series);
            end

            behavior_NWB.nwbdatainterface.set(behavioral_Label,the_position_field_NWB);

            time_epochs(iFile).start_times = behavioral_timestamps(1,1);
            time_epochs(iFile).stop_times  = behavioral_timestamps(end,1);
            time_epochs(iFile).label       = behavioral_Label;
        end

        behavior_NWB.description = 'Behavioral signals';
        nwb.processing.set('behavior', behavior_NWB);

        disp('Behavioral info added..')

    else
        disp('No Yuta behavioral signals present in this directory (*position*)')
        return
    end
end