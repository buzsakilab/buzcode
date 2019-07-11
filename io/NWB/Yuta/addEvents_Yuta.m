function nwb = addEvents_Yuta(xml,nwb)
    %% Add events: nwb2.stimulus_presentation
    eventFiles = dir([xml.folder_path filesep '*.evt']);

    for iFile = 1:length(eventFiles)
        events = LoadEvents(fullfile(eventFiles(iFile).folder,eventFiles(iFile).name)); % LoadEvents is a function from the Buzcode - THIS DOESN'T LOAD ANYTHING FOR 'PulseStim_0V_10021ms_LD0' - maybe because it has a single entrty???
        if ~isempty(events.time) && ~isempty(events.description)
            AnnotationSeries = types.core.AnnotationSeries('data',events.description,'timestamps',events.time);
            nwb.stimulus_presentation.set(events.description{1}, AnnotationSeries);
        end
    end
    disp('Events added..')
end