function nwb = addEpochs_Yuta(xml,nwb)           
            
    %% Add Epochs
    % The epochs are based on separate behavioral files
    % Each separate behavioral file = 1 epoch

    behavioralFiles = dir([xml.folder_path filesep '*Position*']);

    if length(behavioralFiles) == 0
        disp('There are no *Position* files in this folder. No Epochs will be added')
        return
    end


    intervals_epochs = types.core.TimeIntervals();
    id_epochs  = types.core.ElementIdentifiers('data',1:length(behavioralFiles));

    start_time_epochs = zeros(length(behavioralFiles),1); % Start time
    stop_time_epochs  = zeros(length(behavioralFiles),1); % Stop time
    colnames          = cell(length(behavioralFiles),1); % Labels

    for iFile = 1:length(behavioralFiles)

        % The label of the behavior
        behavioral_Label = strsplit(behavioralFiles(iFile).name,'__');
        behavioral_Label = erase(behavioral_Label{2},'.mat');

        position_signals = load(fullfile(behavioralFiles(iFile).folder,behavioralFiles(iFile).name));

        % Some behavioral signal files might have more than one
        % type of signals in them (e.g. twhl_linearized, twhl_norm)
        field_names = fieldnames(position_signals);

        % NO NEED TO KEEP THE DATA - COMMENTING OUT
%                 the_position_field_NWB = types.untyped.Set();
%                 for iField = 1:length(field_names)
%                     the_position_field_NWB.set(field_names{iField}, types.core.SpatialSeries('description', [field_names{iField} ' position signals from ' behavioral_Label ' epoch'], 'data', position_signals.(field_names{iField})(:,2:end)));
%                 end
%                 intervals_epochs.vectordata.set(behavioral_Label,types.core.VectorData('description', ['Position signals from ' behavioral_Label ' epoch'], 'data', the_position_field_NWB));

        start_time_epochs(iFile) = position_signals.(field_names{1})(1,1);
        stop_time_epochs(iFile)  = position_signals.(field_names{1})(end,1);
        colnames{iFile}          = behavioral_Label;

    end

    % Transform doubles to VectorData
    start_time_epochs = types.core.VectorData('description','Starting timestamp of epochs', 'data', start_time_epochs);
    stop_time_epochs  = types.core.VectorData('description','Stopping timestamp of epochs', 'data', stop_time_epochs);

    intervals_epochs.start_time  = start_time_epochs;
    intervals_epochs.stop_time   = stop_time_epochs;
    intervals_epochs.description = 'Experimental epochs';
    intervals_epochs.id          = id_epochs;
    intervals_epochs.colnames    = colnames;

    nwb.intervals_epochs = intervals_epochs;

    disp('Epoch info added..')
end