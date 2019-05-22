function nwb = addEvents(xml,nwb)
    % Add events: nwb2.stimulus_presentation
    % This code was tested on the Buzcode tutorial dataset: 20170505_396um_0um_merge
    % Konstantinos Nasiotis 2019
    
    % Info about IntervalSeries here: https://nwb-schema.readthedocs.io/en/latest/format.html#sec-intervalseries
    
    % All the timestamps enter in a vectorized way (even for 2 dimensional
    % events with start and stop).
    % the start and stop flags are indicated in the .data field by a
    % positive (start) or negative (stop) value
    
    %%
    eventFiles = dir([xml.folder_path filesep '*.events.mat']);
    
    if ~isempty(eventFiles)
        all_events = types.core.ProcessingModule();

        for iFile = 1:length(eventFiles)
            events = load(fullfile(eventFiles(iFile).folder,eventFiles(iFile).name));
            
            name = fields(events); name = name{1}; % Event name (e.g. ripples)
            
            % Get the timestamps
            timestamps = events.(name).timestamps;
       
            if size(timestamps,2) == 2 % Start-stop stimulation - The events need to be saved as IntervalSeries
                IntervalSeries = types.core.IntervalSeries();
                IntervalSeries.timestamps = [timestamps(:,1) ; timestamps(:,2)]; % I linearize the start and stop timestamps: First event start then event stop
                IntervalSeries.data = [ones(size(timestamps,1),1) ; ones(size(timestamps,1),1)* (-1)]; % First event start then event stop: 1 signifies events start, -1 event stop
                
                all_events.nwbdatainterface.set(name,IntervalSeries);
                
                
            else % Stimulation Onset only - The events need to be saved as AnnotationSeries
                AnnotationSeries = types.core.AnnotationSeries('data', repmat({name},length(timestamps),1),'timestamps',timestamps);                
                all_events.nwbdatainterface.set(name,AnnotationSeries);
   
            end
            
            
            
            
% % % % % % % % % % % % %             %% Add the detectorinfo information
% % % % % % % % % % % % %             
% % % % % % % % % % % % %             % This is no properly set. It simply saves the extra info from
% % % % % % % % % % % % %             % the behavior.mat file within the nwb file.
% % % % % % % % % % % % %             % Lab members should help in which information should be saved.
% % % % % % % % % % % % %             
% % % % % % % % % % % % %             detectorinfo = types.core.DynamicTable;
% % % % % % % % % % % % %             
% % % % % % % % % % % % %             detector_fields = fields(events.(name).detectorinfo);
% % % % % % % % % % % % %             for iField = 1:length(detector_fields)
% % % % % % % % % % % % %                    
% % % % % % % % % % % % %                     detectorinfo.vectordata.set(detector_fields{iField}, types.core.VectorData('description','Name of the buzcode function that created the events', 'data', events.(name).detectorinfo.(detector_fields{iField})));
% % % % % % % % % % % % % %                     
% % % % % % % % % % % % % %                     IntervalSeries.timestamps = [events.(name).detectorinfo.(detector_fields{iField})(:,1) ; events.(name).detectorinfo.(detector_fields{iField})(:,2)]; % First event start then event stop
% % % % % % % % % % % % % %                     IntervalSeries.data = [ones(size(events.(name).detectorinfo.(detector_fields{iField}),1),1) ; ones(size(events.(name).detectorinfo.(detector_fields{iField}),1),1)* (-1)]; % First event start then event stop - 1 signifies events start, -1 event stop
% % % % % % % % % % % % % %                     detectorinfo.set(detector_fields{iField}, IntervalSeries);
% % % % % % % % % % % % % %                 else
% % % % % % % % % % % % % %                 	detectorinfo.set(detector_fields{iField}, events.(name).detectorinfo.(detector_fields{iField}));
% % % % % % % % % % % % % %                 end
% % % % % % % % % % % % %             end
% % % % % % % % % % % % %             detectorinfo.colnames = detector_fields;
% % % % % % % % % % % % %             detectorinfo.description = ['Structure fields from ' name ' .behavior.mat file'];
% % % % % % % % % % % % %             detectorinfo.id = types.core.ElementIdentifiers('data', 1);
% % % % % % % % % % % % %             
% % % % % % % % % % % % %             all_events.nwbdatainterface.set([name '_detectorinfo'],detectorinfo);
            

            all_events.description = 'Events created from analysis functions';
            nwb.processing.set('events', all_events);
            
            
        end
        disp('Events added..')
    else
        disp('No *.events.mat files found')
    end
end




