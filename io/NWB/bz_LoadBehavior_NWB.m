function behavior = bz_LoadBehavior_NWB(varargin)

%% Load Behavior from NWB files following the Buzcode format
% behavior = bz_LoadBehavior_NWB;                                                                % Checks the current directory for a NWB file. Gives a list with the available behavior fields in the NWB file for the user to choose which one to load
% behavior = bz_LoadBehavior_NWB('nwb_file', nwb_file);                                          % Checks the specified NWB file. Gives a list with the available behavior fields in the NWB file for the user to choose which one to load
% behavior = bz_LoadBehavior_NWB('nwb_file', nwb_file, 'EightMazePosition_norm_spatial_series'); % Loads the specific behavior from the specified NWB file
% behavior = bz_LoadBehavior_NWB('EightMazePosition_norm_spatial_series');                       % Loads the specific behavior from a NWB file in the current directory


% The behaviors are loaded straight from the NWB file.

% This code creates a file with the behavior selected the first time it's used.
% This is done in order to speed up the process of calling it later 
% (.map takes a long time to compute if the trials are long).


% This code is based on the bz_LoadBehavior
% Konstantinos Nasiotis 2019


%% TODO

disp('1. Check if the BAD trials should be included in the computation of the .map field')




%%
    p = inputParser;
    addParameter(p,'nwb_file','',@isstr);
    
    if size(varargin,2)>2
        behaviorName = varargin{3};
        varargin = {varargin{1,1:2}};
        parse(p,varargin{:})
        nwb_file = p.Results.nwb_file;
    elseif size(varargin,2)==1
        nwb_file   = [];
        behaviorName = varargin{1};
    elseif size(varargin,2)==2
        parse(p,varargin{:})
        nwb_file     = p.Results.nwb_file;
        behaviorName = [];
    elseif size(varargin,2)==0
        nwb_file     = [];
        behaviorName = [];
    end


    %% let's check that there is an appropriate NWB file
    if isempty(nwb_file)
       %disp('No nwb_file given, so we look for a *nwb file in the current directory')
       d = dir('*nwb');
       if length(d) > 1 % we assume one .nwb file or this should break
           error('there is more than one .nwb file in this directory?');
       elseif length(d) == 0
           d = dir('*nwb');
           if isempty(d)
               error('could not find an NWB file..')
           end
       end
       nwb_file = fullfile(d.folder, d.name);
    end
    

    nwb2 = nwbRead(nwb_file);

    % Check if behavior fields exists in the dataset
    try
        allBehaviorKeys = keys(nwb2.processing.get('behavior').nwbdatainterface);
        
        behavior_exist_here = ~isempty(allBehaviorKeys);
        if ~behavior_exist_here
            disp('No behavior in this .nwb file')
            return
            
        else
            disp(' ')
            disp('The following behavior types are present in this dataset')
            disp('------------------------------------------------')
            for iBehavior = 1:length(allBehaviorKeys)
                disp(allBehaviorKeys{iBehavior})
            end
            disp(' ')
        end
    catch
        disp('No behavior in this .nwb file')
        return
    end

    %% Check if the behavior file has already been computed (it takes some time to fill the .mapping field, especially when the trials are long)
    %  If the behavior file is already computed, just load the file and return
    
    % FIRST CHECK - CHECK IF THE behaviorName REQUESTED EXISTS
    
    [the_path, name, ext] = fileparts(nwb_file);
    
    if exist(fullfile(the_path, [name '.' behaviorName '.behavior.mat']) , 'file') == 2
        load(fullfile(the_path, [name '.' behaviorName '.behavior.mat']))
        disp('Just loading the already computed behavior struct')
        return
    else
        disp('A file with the specified Behavior hasnt already been created. Select the behavior from the list')
    end
    
    %% If a the behavior file doesn't exist, create a new one
    
    % Check if a specific behavior was called to be loaded. If not, display a
    % pop-up list with the available behaviors for selection
    [iBehavior, ~] = listdlg('PromptString','Which behavior type would you like to load?',...
                             'ListString',allBehaviorKeys,'SelectionMode','single');
  
    % Fill the behavior 
    behavior = struct; % Initialize

    
    % Select the key that is within the Behavior selection
    allBehaviorSUBKeys = keys(nwb2.processing.get('behavior').nwbdatainterface.get(allBehaviorKeys(iBehavior)).spatialseries);    
    if length(allBehaviorSUBKeys)>1
        [iSubKeyBehavior, ~] = listdlg('PromptString','Multiple Behavior channel-selections within the selected behavior',...
                                     'ListString',allBehaviorSUBKeys,'SelectionMode','single');
    else
        iSubKeyBehavior = 1;
    end
        
    %% Check if the behavior file has already been computed (it takes some time to fill the .mapping field, especially when the trials are long)
    %  If the behavior file is already computed, just load the file and return

    % SECOND CHECK - HERE A SUBKEY HAS ALREADY BEEN SET FROM THE SELECTION
    
    if exist(fullfile(the_path, [name '.' allBehaviorSUBKeys{iSubKeyBehavior} '.behavior.mat']) , 'file') == 2
        load(fullfile(the_path, [name '.' allBehaviorSUBKeys{iSubKeyBehavior} '.behavior.mat']))
        disp('Just loading the already computed behavior struct')
        return
    end
    
    
    %% Get the behavior info
    disp('computing')
    
    selected_behavior = nwb2.processing.get('behavior').nwbdatainterface.get(allBehaviorKeys(iBehavior)).spatialseries.get(allBehaviorSUBKeys{iSubKeyBehavior}); % This is for easier reference through the code
    
    behavior.samplingRate = 1000 ./ nanmean(diff(selected_behavior.timestamps.load))./1000;
    behavior.units        = selected_behavior.data_unit;

    Info = allBehaviorKeys{iBehavior};
%     behavior.description  = selected_behavior.description;

    % Check if timestamps and data have values. If not something weird
    % is going on
    if ~isempty(selected_behavior.timestamps)
        behavior.timestamps = selected_behavior.timestamps.load;
    else
        behavior.timestamps     = [];
        warning(['Behavior: ' Info ' --- Timestamps are empty: weird'])
    end
    if ~isempty(selected_behavior.data)
        
        % Specifically for position signals, seperate to x and y
        if ~isempty(strfind(allBehaviorKeys(iBehavior),'position'))
            data = selected_behavior.data.load;
            datasubstruct.x = data(1,:)';
            if size(data,1)>1
                datasubstruct.y = data(2,:)';
            end
            if size(data,1)>2
                datasubstruct.z = data(3,:)';
            end
        else
            datasubstruct.data = selected_behavior.data.load;
            
        end
        
        
    else
        datasubstruct.data = [];
        warning(['Behavior: ' Info ' --- Data is empty: weird'])
    end


    % position: .x, .y, and .z
    % units: millimeters, centimeters, meters[default]
    % orientation: .x, .y, .z, and .w
    % rotationType: euler or quaternion[default]
    % pupil: .x, .y, .diameter


%       datasubstruct.reference_frame     = nwb2.acquisition.get(all_raw_keys(allBehaviorKeys(iBehavior))).reference_frame; This seems to be present only on the position_sensor channels
   
    % Conside removing these
    datasubstruct.starting_time_unit  = selected_behavior.starting_time_unit;
    datasubstruct.timestamps_interval = selected_behavior.timestamps_interval;
    datasubstruct.comments            = selected_behavior.comments;
    datasubstruct.control             = selected_behavior.control;
    datasubstruct.control_description = selected_behavior.control_description;
    datasubstruct.data_resolution     = selected_behavior.data_resolution;
    datasubstruct.starting_time       = selected_behavior.starting_time;
    datasubstruct.help                = selected_behavior.help;

    
    
    %% Exclude special characters from the behavior name
    BehaviorLabel = allBehaviorKeys{iBehavior};    
    
    clean_string = regexprep(BehaviorLabel,'[^a-zA-Z0-9_]','');
    if ~strcmp(clean_string, BehaviorLabel)
        disp(['The variable name (' BehaviorLabel ') of the Behavior was changed to exclude special characters'])
        BehaviorLabel = clean_string;
    end
    
    %% If this is a position signal, rename to position to have compatibility with the Buzsaki functions
    if ~isempty(strfind(allBehaviorKeys(iBehavior),'position'))
        BehaviorLabel = 'position';
    end
    
    
    behavior.(BehaviorLabel)    = datasubstruct;
    
    % BehaviorInfo
    behaviorinfo.description       = selected_behavior.description;
    behaviorinfo.acquisitionsystem = 'Fill Me';
    behaviorinfo.substructnames    = {BehaviorLabel};
    behavior.behaviorinfo          = behaviorinfo;


    
    %% Fill the events substructure
    
    nTrials          = length(nwb2.intervals_trials.start_time.data.load);
    uniqueConditions = unique(nwb2.intervals_trials.vectordata.get('condition').data);
    
    
    all_conditions = nwb2.intervals_trials.vectordata.get('condition').data;
    for iCondition = 1:length(all_conditions)
        
        events.trialConditions(iCondition) = find(strcmp(all_conditions{iCondition}, uniqueConditions));%1x221 double     which unique condition each trialCondition belongs to - INDEX

    end
    events.trialIntervals  = [nwb2.intervals_trials.start_time.data.load nwb2.intervals_trials.stop_time.data.load];       % 221x2 double   start-end of trial in seconds
    events.conditionType   = unique(nwb2.intervals_trials.vectordata.get('condition').data)';     % 1x10 cell   (central, wheel ...)
    
%     nwb2.intervals_trials.id.data.load
%     nwb2.intervals_trials.colnames
%     nwb2.intervals_trials.vectordata
%     nwb2.intervals_trials.vectordata.get('both_visit').data.load
%     nwb2.intervals_trials.vectordata.get('condition').data
%     nwb2.intervals_trials.vectordata.get('error_run').data.load
%     nwb2.intervals_trials.vectordata.get('stim_run').data.load
%     nwb2.intervals_trials.vectorindex
    
    trialTimestampBounds = [nwb2.intervals_trials.start_time.data.load nwb2.intervals_trials.stop_time.data.load];
    position_timestamps =  selected_behavior.timestamps.load;

    for iTrial = 1:nTrials
        
        %% Load only the samples from each trial
        % Find the indices of the timestamps that are selected
        
        [~, iPositionTimestamps] = histc(trialTimestampBounds(iTrial,:), position_timestamps);
        
        if iPositionTimestamps(1)~=0 && iPositionTimestamps(2)==0
            iPositionTimestamps(2) = length(position_timestamps);
        end

        try
            position = selected_behavior.data.load([1, iPositionTimestamps(1)], [Inf, iPositionTimestamps(2)]);
        catch
            error ('Error loading selected behavior file. Are there datapoints within the trials selection time-period?')
        end
        
        % The boundaries struct will hold the position boundaries of each
        % TYPE of CONDITION, for all the trials on each condition Type
        if ~exist('boundaries')
            if size(position,1)==1
                boundaries = repmat(struct('x_min',0,'x_max',0),length(uniqueConditions),1);
            elseif size(position,1)==2
                boundaries = repmat(struct('x_min',0,'x_max',0,'y_min',0,'y_max',0),length(uniqueConditions),1);
            elseif size(position,1)==3
                boundaries = repmat(struct('x_min',0,'x_max',0,'y_min',0,'y_max',0,'z_min',0,'z_max',0),length(uniqueConditions),1);    
            end
        end

        
%         clr = rand(1,3);
%         plot(position(1,:),position(2,:),'.','color',clr);
%         drawnow
                
        % Get the index of the condition that the trials belong to - this will
        % be used for defining the spatial boundaries of the condition
        iConditionOfTrial = find(strcmp(uniqueConditions, nwb2.intervals_trials.vectordata.get('condition').data(iTrial)))';

        
        trials{1,iTrial}.x = position(1,:)';% 608 x 1
        if size(position,1)>1
            trials{1,iTrial}.y = position(2,:)';% 608 x 1
        end
        if size(position,1)>2
            trials{1,iTrial}.z = position(3,:)';% 608 x 1
        end
        
        trials{1,iTrial}.timestamps = position_timestamps(iPositionTimestamps(1):iPositionTimestamps(2));% 608 x 1
        trials{1,iTrial}.direction  = 'FILL ME'; % 'clockwise' 'counterclockwise'
        trials{1,iTrial}.type       = nwb2.intervals_trials.vectordata.get('condition').data{iTrial}; % 'central alternation'
    %     trials{1,iTrial}.orientation    = % 1x1 struct
    %     trials{1,iTrial}.errorPerMarker = 608 x 1

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Map the trials to a 1x201 vector     
    bins = 200;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % normalize positions to template
    c=1;

    uniqueConditions = unique((events.trialConditions));

    for iCondition = 1:length(uniqueConditions)

        map{iCondition}=[];
        t_conc=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MAYBE SELECT ONLY THE TRIALS THAT WEREN'T BAD HERE
%         iTrialsInCondition = find(events.trialConditions == uniqueConditions(iCondition) & ~nwb2.intervals_trials.vectordata.get('error_run').data.load');
        iTrialsInCondition = find(events.trialConditions == uniqueConditions(iCondition));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        

        for iTrial = 1:length(iTrialsInCondition)
            selected_trial = trials{iTrialsInCondition(iTrial)}; % I set this for faster typing

            if size(position,1)==1
                t_conc = [selected_trial.timestamps, selected_trial.x, 20*(selected_trial.timestamps - selected_trial.timestamps(1))]; % Check what this 20 is
            elseif size(position,1)==2
                t_conc = [selected_trial.timestamps, selected_trial.x selected_trial.y, 20*(selected_trial.timestamps - selected_trial.timestamps(1))]; % Check what this 20 is
            elseif size(position,1)==3
                t_conc = [selected_trial.timestamps, selected_trial.x selected_trial.y selected_trial.z, 20*(selected_trial.timestamps - selected_trial.timestamps(1))]; % Check what this 20 is
            end
            
% % % % % %     templateVector = 1:nBins;
% % % % % %     target = size(position,2);
% % % % % %     n_ent = length(templateVector);
% % % % % %     trials{1,iTrial}.mapping = round(interp1( 1:n_ent, templateVector, linspace(1, n_ent, target) ))';

            disp('the computation below doesnt make the usage of this function practical for on the fly retrieval of the Behavior')
            if length(t_conc)>=bins
                while length(t_conc)>bins+1      
                    
                    di = pdist(t_conc);
                    s = squareform(di);
                    s(find(eye(size(s))))=nan;
                    [a b] = min(s(:));
                    [coords blah] = find(s==a);
                    t_conc(coords(1),:) = (t_conc(coords(1),:)+t_conc(coords(2),:))./2;
                    t_conc(coords(2),:) = [];
                end
                t_conc_all(iTrial,:,:) = t_conc;
            end
        end
        if length(iTrialsInCondition)>0
            map{iCondition} = squeeze(median(t_conc_all(:,:,:),1));
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assign to the right behavior field
        if size(position,1)==1 % 1D
            events.map{iCondition}.x = map{iCondition}(:,2);
        elseif size(position,1)==2 % 2D   
            events.map{iCondition}.x = map{iCondition}(:,2);
            events.map{iCondition}.y = map{iCondition}(:,3);
        elseif size(position,1)==3 % 3D
            events.map{iCondition}.x = map{iCondition}(:,2);
            events.map{iCondition}.y = map{iCondition}(:,3);
            events.map{iCondition}.z = map{iCondition}(:,4);   
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clear t_conc_all

        disp('finding mapping...')
        for iTrial =1:length(iTrialsInCondition)  % all trial types (rotations)
            selected_trial = trials{iTrialsInCondition(iTrial)}; % I set this for faster typing
            for iPoint = 1:length(selected_trial.timestamps)
                
                if size(position,1)==1 % 1D
                    [a b] = min(nansum(abs([selected_trial.timestamps(iPoint)-map{iCondition}(:,1),...   % Timestamp
                                            selected_trial.x(iPoint)-map{iCondition}(:,2),...   % X_COORDINATE
                                           (selected_trial.timestamps(iPoint)-selected_trial.timestamps(1))*50-map{iCondition}(:,1),...  % penalty for time differences
                                            (40*(iPoint./length(selected_trial.timestamps)*length(map{iCondition}) - (1:length(map{iCondition})))')])'));     % penalty for order differences
                    
                elseif size(position,1)==2 % 2D
                    [a b] = min(nansum(abs([selected_trial.timestamps(iPoint)-map{iCondition}(:,1),...   % Timestamp
                                            selected_trial.x(iPoint)-map{iCondition}(:,2),...   % X_COORDINATE
                                            selected_trial.y(iPoint)-map{iCondition}(:,3),...   % Y_COORDINATE 
                                           (selected_trial.timestamps(iPoint)-selected_trial.timestamps(1))*50-map{iCondition}(:,1),...  % penalty for time differences
                                            (40*(iPoint./length(selected_trial.timestamps)*length(map{iCondition}) - (1:length(map{iCondition})))')])'));     % penalty for order differences
                    
                elseif size(position,1)==3 % 3D
                    [a b] = min(nansum(abs([selected_trial.timestamps(iPoint)-map{iCondition}(:,1),...   % Timestamp
                                            selected_trial.x(iPoint)-map{iCondition}(:,2),...   % X_COORDINATE
                                            selected_trial.y(iPoint)-map{iCondition}(:,3),...   % Y_COORDINATE 
                                            selected_trial.z(iPoint)-map{iCondition}(:,4),...   % Z_COORDINATE
                                           (selected_trial.timestamps(iPoint)-selected_trial.timestamps(1))*50-map{iCondition}(:,1),...  % penalty for time differences
                                            (40*(iPoint./length(selected_trial.timestamps)*length(map{iCondition}) - (1:length(map{iCondition})))')])'));     % penalty for order differences
                end
                                    
                mapping{iCondition}{iTrial}(iPoint,:) = [map{iCondition}(b,1:end) b selected_trial.timestamps(iPoint)];
            
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Assign to the right behavior field
            trials{iTrialsInCondition(iTrial)}.mapping = mapping{iCondition}{iTrial}(:,end-1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    
    events.trials   = trials;
    behavior.events = events;
    
    %% Check that the behavior structure meets buzcode standards
    [isBehavior] = bz_isBehavior(behavior);
    switch isBehavior
        case false
            warning('Your behavior structure does not meet buzcode standards. Sad.')
    end
    
    %% Save the behavior files if they are not already saved
    if exist(fullfile(the_path, [name '.' allBehaviorSUBKeys{iSubKeyBehavior} '.behavior.mat']) ,'file') ~=2
        save(fullfile(the_path, [name '.' allBehaviorSUBKeys{iSubKeyBehavior} '.behavior.mat']), 'behavior')
    end
    
end