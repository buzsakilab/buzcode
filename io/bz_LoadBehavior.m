function [ behavior,filename ] = bz_LoadBehavior( basePath,behaviorName )
%[ behavior,filename ] = bz_LoadBehavior( basePath,behaviorName ) function for
%loading behavior.mat files.
%
%
%DLevenstein 2017

% Added NWB support: Konstantinos Nasiotis 2019
%%

if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);

if ~exist('behaviorName','var')
    allBehaviorFiles = dir(fullfile(basePath,[baseName,'.','*','.behavior.mat']));
    
    if ~isempty(allBehaviorFiles)
        [s,v] = listdlg('PromptString','Which behavior.mat would you like to load?',...
                     'ListString',{allBehaviorFiles.name},'SelectionMode','single');
        if isempty(s)
            behavior = []; filename = [];
            return
        end
        filename = fullfile(basePath,allBehaviorFiles(s).name);
    else
        check_for_NWB = 1;
        filename = [];
    end
        
else
    filename = fullfile(basePath,[baseName,'.',behaviorName,'.behavior.mat']);
end


if exist(filename,'file')
    behaviorstruct = load(filename);
    check_for_NWB = 0;
else
    check_for_NWB = 1;
end


if check_for_NWB
    if ~exist('behaviorName','var')
        behaviorName = [];        
    end
    behaviorstruct = loadBehaviorNWB(basePath, behaviorName);
    if isempty(behaviorstruct)
       behavior = []; filename = [];
       return
    end
end



varsInFile = fieldnames(behaviorstruct);

if numel(varsInFile)==1
    behavior = behaviorstruct.(varsInFile{1});
else
    warning('Your .behavior.mat has multiple variables/structures in it... wtf.')
    behavior = behaviorstruct;
end

%Check that the behavior structure meets buzcode standards
[isBehavior] = bz_isBehavior(behavior);
switch isBehavior
    case false
        warning('Your behavior structure does not meet buzcode standards. Sad.')
end





%% Function that loads NWB fields to a Behavior structure
function behaviorstruct = loadBehaviorNWB(basePath, behaviorName)
    
    % Check for NWB file since no .behavior.mat files were found
    d = dir(fullfile(basePath, '*.nwb'));
    if length(d) > 1 % If more than one .nwb files exist in the directory, select which one to load from
        warning('there is more than one .nwb file in this directory');
        [iNWBFile, ~] = listdlg('PromptString','Which NWB file would you like to load?',...
                                 'ListString',{d.name},'SelectionMode','single');
        d = d(iNWBFile);
    elseif length(d) == 0
        d = dir('*nwb');
        if isempty(d)
            disp('No .behavior.mat or NWB files present in this directory')
            behaviorstruct = [];
            return
        end
    end
    nwb_file = fullfile(d.folder, d.name); 


    
    if exist(nwb_file,'file') == 2
        nwb2 = nwbRead(nwb_file);
        
        % Check if the behavior field exists in the dataset
        all_processing_keys   = keys(nwb2.processing);
        behavior_FIELD_exists = sum(ismember(all_processing_keys, 'behavior'))==1;
        
        if behavior_FIELD_exists
            allBehaviorKeys = keys(nwb2.processing.get('behavior').nwbdatainterface);
            
            % Remove the trials_ keys since they are just the trials sub-field
            % within the events field of the behavior.mat files
            ii = [];
            for iKey = 1:length(allBehaviorKeys)
                if isempty(strfind(allBehaviorKeys{iKey}, 'trials_'))
                    ii = [ii,iKey];
                end
            end
            allBehaviorKeys = {allBehaviorKeys{ii}};
            
        
            behavior_exist_here = ~isempty(allBehaviorKeys);
            if ~behavior_exist_here
                disp('No entry is filled in the Behavior field')
                behaviorstruct = [];
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
            
            
            % Check if a specific behavior was called to be loaded. If not, display a
            % pop-up list with the available behaviors for selection
            
            % NOTE: THE PREFIX TRIALS_ WILL BE IGNORED SINCE IT SIMPLY
            % IMPLIES THE TRIALS THAT ARE WITHIN THE BEHAVIOR.MAT FILE
            
            if isempty(behaviorName)
                [iBehavior, ~] = listdlg('PromptString','Which behavior type would you like to load?',...
                                         'ListString',allBehaviorKeys,'SelectionMode','single');
            else
                if sum(ismember(allBehaviorKeys, behaviorName))==1 % If the behavior exists, find its index
                    iBehavior = find(ismember(allBehaviorKeys, behaviorName));
                else
                    disp('The requested behavior doesnt exist in this NWB file')
                    behaviorstruct = [];
                    return
                end
            end
            
            %% Check if the Behavior fields are filled from a file that meets standard Buzcode standards (*.behavior.mat)
            %  This flag is added from the addBehavior function
            
            if isempty(strfind(nwb2.processing.get('behavior').description, 'behavior.mat file'))
                error('The behavior in this nwb file was not added from a file that meets the standard Buzcode format')
            end
            
            %% If the check for standard format passed, get the behavior
            % fields
            
            % Fill the behavior 
            behaviorstruct = struct; % Initialize

            % Select the key that is within the Behavior selection
            allBehaviorfields = keys(nwb2.processing.get('behavior').nwbdatainterface.get(allBehaviorKeys(iBehavior)).spatialseries);    
            
            %% Fill the Behavior Structure - This is added to the NWB file from addBehavior
            
            for iField = 1:length(allBehaviorfields)
            
                if strcmp('behaviorinfo',allBehaviorfields{iField})
                    selected_behavior = nwb2.processing.get('behavior').nwbdatainterface.get(allBehaviorKeys(iBehavior)).spatialseries.get(allBehaviorfields{iField}); % This is for easier reference through the code
                    behaviorstruct.behavior.samplingRate = 1000 ./ nanmean(diff(selected_behavior.timestamps.load))./1000;
                    behaviorstruct.behavior.units        = selected_behavior.data_unit;
                    behaviorstruct.behavior.rotationType = '';
                    
                    behaviorstruct.behavior.timestamps = selected_behavior.timestamps.load;
                    
                    behaviorstruct.behavior.behaviorinfo.description       = selected_behavior.description;
                    behaviorstruct.behavior.behaviorinfo.substructnames    = {allBehaviorfields{~ismember(allBehaviorfields, 'behaviorinfo')}};
                    behaviorstruct.behavior.behaviorinfo.errorPerMarker    = selected_behavior.data.load;
                    behaviorstruct.behavior.behaviorinfo.acquisitionsystem = selected_behavior.control_description.load;
                    
                else
                    selected_behavior = nwb2.processing.get('behavior').nwbdatainterface.get(allBehaviorKeys(iBehavior)).spatialseries.get(allBehaviorfields{iField}); % This is for easier reference through the code

                    axis = {'x','y','z','w'}; % I assume that this is the order that the axis are filled
                    
                    nAxis = selected_behavior.data.dims(1); % 4 x 263958 for orientation
                    
                    for iAxis = 1:nAxis
                        behaviorstruct.behavior.(allBehaviorfields{iField}).(axis{iAxis}) = selected_behavior.data.load([1,iAxis],[Inf,iAxis]);
                    end
                end
            end
            
            %% Fill the Events - Trials subStructure - This is added to the NWB file from addTrials
            trials_object = nwb2.processing.get('behavior').nwbdatainterface.get(['trials_' allBehaviorKeys{iBehavior}]); % Time intervals object
            
            % SECOND CHECK THAT THE TRIALS FIELD HAS BEEN FILLED FROM A
            % BEHAVIOR.MAT FILE
            if isempty(strfind(trials_object.description, 'behavior.mat file'))
                error('The trials in this nwb file was not added from a file that meets the standard Buzcode format')
            end
            
            nTrials = trials_object.start_time.data.dims;
            
            for iTrial = 1:nTrials
                for iColname = 1:length(trials_object.colnames)
                    
                    singleTrialData = trials_object.vectordata.get('trials').data.load([1,iTrial,iColname],[Inf,iTrial,iColname]);
                    singleTrialData = singleTrialData(~isnan(singleTrialData)); % I concatenated with NaNs in the addTrials function in order to use a matrix and not a cell array
                    
                    if isempty(strfind(trials_object.colnames{iColname},'orientation'))
                        behaviorstruct.behavior.events.trials{iTrial}.(trials_object.colnames{iColname}) = singleTrialData;
                    else
                        behaviorstruct.behavior.events.trials{iTrial}.orientation.(trials_object.colnames{iColname}(end)) = singleTrialData;
                    end
                end
                behaviorstruct.behavior.events.trials{iTrial}.direction = trials_object.vectordata.get('direction').data.load(iTrial);
                behaviorstruct.behavior.events.trials{iTrial}.type      = trials_object.vectordata.get('type').data.load(iTrial);
            end
            
            
            for iMap = 1: trials_object.vectordata.get('map').data.dims(2) % Number of conditions
                axis = {'x','y','z'};
                for iAxis = 1:length(axis)
                    behaviorstruct.behavior.events.map{iMap}.(axis{iAxis}) = trials_object.vectordata.get('map').data.load([1,iMap,iAxis], [Inf,iMap,iAxis]);
                end
            end
            
            behaviorstruct.behavior.events.trialIntervals = [trials_object.start_time trials_object.stop_time];
            behaviorstruct.behavior.events.conditionType  = trials_object.vectordata.get('conditionType').data.load;
            
            %% Optional
            % reorder for easy comparison with the tutorial example
            C = {'position','timestamps','samplingRate','units','orientation','rotationType','events','behaviorinfo'};
            behaviorstruct.behavior = orderfields(behaviorstruct.behavior,C);
            
            
        else
            disp('No behavior in this NWB file')
            behaviorstruct = [];
            return
        end
            
    else
        disp('No .behavior.mat or NWB files were found')
        behaviorstruct = [];
        return
    end
    

end



end