function nwb = addTrials(xml,nwb)
            
    % THIS SECTION ADDS THE TRIAL INFO FROM THE .BEHAVIOR.MAT FILES TO THE NWB FILE - COMBINE IT WITH addBehavior in order to use bz_LoadBehavior
    % The trials part of Each file is stored in: nwb.processing.get('behavior').nwbdatainterface.get('trials_BEHAVIORNAME'))
    % This is based on the Buzcode tutorial Behavior file: 20170505_396um_0um_merge.track.behavior.mat
    % and the Buzcode wiki: https://github.com/buzsakilab/buzcode/wiki/Data-Formatting-Standards#behavior
    
    % Konstantinos Nasiotis 2019  
    %% Fill the fields
    behavioralFiles = dir([xml.folder_path filesep '*behavior.mat']);

    if ~isempty(behavioralFiles)

        intervals_trials = types.core.TimeIntervals();

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

            all_start_times = types.core.VectorData('description','Starting timepoint of Each Trial','data',behavior.events.trialIntervals(:,1));
            all_stop_times  = types.core.VectorData('description','Ending timepoint of Each Trial','data',behavior.events.trialIntervals(:,2));

            

            %% Add optional events field
            if isfield(behavior, 'events')

                trials = behavior.events.trials;
                nTrials = length(trials);

                % Find the maximum length of a trial
                % All trials will be concatenated on a single matrix
                % Smaller trials than the maximum will be filled with
                % NaNs

                maxLength = 0;
                for iTrial = 1:nTrials
                    if maxLength<length(trials{iTrial}.x)
                        maxLength = length(trials{iTrial}.x);
                    end
                end

                % Check what fields exist:
                all_trial_fields = fieldnames(trials{iTrial});
                presentFields = all_trial_fields(ismember(all_trial_fields, {'x', 'y', 'z', 'errorPerMarker', 'mapping','orientation', 'timestamps', 'direction', 'type'})');

                % Check orientation fields
                for iField = 1:length(presentFields)
                    % The orientation fields has subfields (other channels)
                    if strcmp(presentFields{iField}, 'orientation')
                        all_Orientation_fields = fieldnames(trials{iTrial}.(presentFields{iField}));
                    end
                end

                isOrientationPresent = ismember('orientation', presentFields);
                isDirectionPresent   = ismember('direction', presentFields);
                isTypePresent        = ismember('type', presentFields);

                oneMatrixToRuleThemAll = zeros(maxLength, nTrials, length(presentFields) + isOrientationPresent* (length(all_Orientation_fields) - 1) - isDirectionPresent - isTypePresent);
                labels    = cell(size(oneMatrixToRuleThemAll, 3), 1); % the labels of the vectors. It will imply the order they are saved on the matrix
                direction = cell(nTrials,1);
                type      = cell(nTrials,1);

                % Fill the matrix with all vectors
                for iTrial = 1:nTrials
                    ii = 1;
                    for iVector = 1:length(presentFields)
                        if strcmp(presentFields{iVector}, 'direction')
                            direction{iTrial} = trials{iTrial}.direction;

                        elseif strcmp(presentFields{iVector}, 'type')
                            type{iTrial} = trials{iTrial}.type;

                        elseif strcmp(presentFields{iVector}, 'orientation')
                            for iOrientation = 1:length(all_Orientation_fields)
                                labels{ii} = ['orientation_' all_Orientation_fields{iOrientation}];

                                % I CONCATENATE WITH NANS
                                oneMatrixToRuleThemAll(:,iTrial,ii) = [trials{iTrial}.orientation.(all_Orientation_fields{iOrientation}); zeros(maxLength - length(trials{iTrial}.orientation.(all_Orientation_fields{iOrientation})),1)*NaN];
                                ii = ii + 1;
                            end
                        else
                            labels{ii} = presentFields{iVector};
                            oneMatrixToRuleThemAll(:,iTrial,ii) = [trials{iTrial}.(presentFields{iVector}); zeros(maxLength - length(trials{iTrial}.(presentFields{iVector})),1)*NaN];
                            ii = ii + 1;
                        end
                    end
                end

                % To sum up, the oneVectorToRuleThemAll holds all info vectors from the trials. The labels variable holds the labels for each one
                % Add in one MATRIX the .map info as well

                all_map_fields = fieldnames(behavior.events.map{1});
                map_matrix = zeros(length(behavior.events.map{1}.(all_map_fields{1})),length(behavior.events.map), length(all_map_fields));
                for iCondition = 1:length(behavior.events.map)
                    for iMapField = 1:length(all_map_fields)
                        map_matrix(:,iCondition,iMapField) = behavior.events.map{iCondition}.(all_map_fields{iMapField});
                    end
                end

                conditionType = behavior.events.conditionType;

                %% Dump the vectors from all trials to the nwb
                intervals_trials.vectordata.set('trials', types.core.VectorData('description', ['all trials" data from ' behavioral_Label '.behavior.mat file - Shorter trials were concatenated with NaNs - nSamples x nTrials x nChannels'], 'data', oneMatrixToRuleThemAll));
                intervals_trials.vectordata.set('map', types.core.VectorData('description', ['map info from all condition for ' behavioral_Label '.behavior.mat file - nSamples x nConditions x nChannels'], 'data', map_matrix));
                intervals_trials.vectordata.set('conditionType', types.core.VectorData('description', ['conditionType for map field in ' behavioral_Label], 'data', conditionType));
                intervals_trials.vectordata.set('direction', types.core.VectorData('description', ['direction for trials in ' behavioral_Label], 'data', direction));
                intervals_trials.vectordata.set('type', types.core.VectorData('description', ['Type of trials in ' behavioral_Label], 'data', type));

                intervals_trials.colnames = labels;
                intervals_trials.id       = types.core.ElementIdentifiers('data', 1:nTrials);
                intervals_trials.start_time = all_start_times;
                intervals_trials.stop_time = all_stop_times;                    

                intervals_trials.description = ['Trial info from ' behavioral_Label '.behavior.mat file']; % use the description for differentiation with the other condition (loading from several .mat files)

            end

            % Check if the behavior field is already created           
            if isempty(keys(nwb.processing))
                behavior_NWB = types.core.ProcessingModule('description', 'intermediate data from extracellular electrophysiology recordings, e.g., LFP');
                nwb.processing.set('behavior', behavior_NWB);
            else
                if ~ismember(keys(nwb.processing),'behavior')
                    behavior_NWB = types.core.ProcessingModule('description', 'intermediate data from extracellular electrophysiology recordings, e.g., LFP');
                    nwb.processing.set('behavior', behavior_NWB);
                end
            end

            nwb.processing.get('behavior').nwbdatainterface.set(['trials_' behavioral_Label], intervals_trials);

        end

        disp('Trial info added..')

    else
        disp('No *behavior.mat file present in the selected directory..')
    end


end