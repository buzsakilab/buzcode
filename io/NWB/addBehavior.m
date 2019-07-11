function nwb = addBehavior(xml,nwb)
    %% This function adds Behavior information on an nwb file - COMBINE IT WITH addTrials in order to use bz_LoadBehavior
    %  It is based on the Buzcode tutorial Behavior file: 20170505_396um_0um_merge.track.behavior.mat
    %  and the Buzcode wiki: https://github.com/buzsakilab/buzcode/wiki/Data-Formatting-Standards#behavior

    % Konstantinos Nasiotis 2019
    %% Add behavioral data: nwb2.processing.get('behavior').nwbdatainterface
    behavior_NWB = types.core.ProcessingModule;

    behavioralFiles = dir([xml.folder_path filesep '*behavior.mat']);

    if length(behavioralFiles) ~= 0

        for iFile = 1:length(behavioralFiles)

            % The label of the behavior
            behavioral_Label = erase(behavioralFiles(iFile).name,'.behavior.mat');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            behavioral_Label = strsplit(behavioral_Label,'.'); % '20170505_396um_0um_merge'    'track'
            behavioral_Label = behavioral_Label{2};            % This section is hardcoded, maybe improve.
                                                               % I assumed here that the standardized way of saving behavior files is: experimentName.BehaviorLabel.behavior.mat
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Load the file
            behaviorstruct = load(behavioralFiles(iFile).name); % The example I have has the variable "behavior" saved
            behavior       = behaviorstruct.behavior;


            
            behavioral_timestamps = behavior.timestamps;


            behavioral_signals_NWB = types.core.Position;
            behavior_field_names = fieldnames(behaviorstruct.behavior)';

            %% First add the fields that contain channels
            iChannelTypeFields = find(ismember(behavior_field_names,{'position','orientation','pupil'}));

            for iField = iChannelTypeFields

                channel_fields      = fieldnames(behavior.(behavior_field_names{iField}));
                behavioral_channels = zeros(length(behavior.(behavior_field_names{iField}).(channel_fields{1})),length(channel_fields));

                for iChannel = 1:length(channel_fields)
                    behavioral_channels(:,iChannel) = behavior.(behavior_field_names{iField}).(channel_fields{iChannel});
                end

                spatial_series = types.core.SpatialSeries('data', behavioral_channels, 'timestamps', behavioral_timestamps,'data_unit', behavior.units);
                behavioral_signals_NWB.spatialseries.set(behavior_field_names{iField}, spatial_series);
            end
            behavior_NWB.nwbdatainterface.set(behavioral_Label,behavioral_signals_NWB);


            %% Add behaviorinfo field - THIS SHOULDN'T BE SPATIALSERIES - HOWEVER IF IT IS NOT, THE IMPORTING FAILS
            behaviorinfoField = find(ismember(behavior_field_names,{'behaviorinfo'}));

            behaviorInfo_fields = fieldnames(behavior.(behavior_field_names{behaviorinfoField}));

            ErrorPerMarkerSignal = behavior.(behavior_field_names{behaviorinfoField}).errorPerMarker;

            spatial_series = types.core.SpatialSeries('data', ErrorPerMarkerSignal, 'timestamps', behavioral_timestamps,'data_unit', behavior.units, ...
                                                      'comments', 'The data field represent the errorPerMarker vector', 'description', behavior.(behavior_field_names{behaviorinfoField}).description, 'control_description', behavior.(behavior_field_names{behaviorinfoField}).acquisitionsystem);
            behavioral_signals_NWB.spatialseries.set(behavior_field_names{behaviorinfoField}, spatial_series);
            behavior_NWB.nwbdatainterface.set(behavioral_Label,behavioral_signals_NWB);

        end

        behavior_NWB.description = ['Behavioral signals from ' behavioral_Label 'behavior.mat file'];
        nwb.processing.set('behavior', behavior_NWB);
        disp('Behavioral info added..')

    else
        disp('No behavior.mat file present in this directory')
        return
    end
end