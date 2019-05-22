function [ events,filename ] = bz_LoadEvents(basePath,eventsName)
%[ events ] = bz_LoadEvents(basePath,eventsName) function for
%loading events.mat files. events.mat files are saved as...
% datasetPath/baseName/baseName.eventsName.events.mat
%
%eventsName can be the name of a events.mat file, If empty, prompts the user
%with a list of available events.mat files in basePath.
%Future update: 'all' (nonfunctional) to load all events.mat files for a given recording. 
%
%DLevenstein 2017

% Added NWB support: Konstantinos Nasiotis 2019
%%
if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);

if ~exist('eventsName','var')
    allEventsFiles = dir(fullfile(basePath,[baseName,'.','*','.events.mat']));
    
    if ~isempty(allEventsFiles)
        [s,v] = listdlg('PromptString','Which events.mat would you like to load?',...
                     'ListString',{allEventsFiles.name},'SelectionMode','single');
        if isempty(s)
            events = []; filename = [];
            return
        end
        filename = fullfile(basePath,allEventsFiles(s).name);
    else
        check_for_NWB = 1;
        filename = [];
    end
else
    filename = fullfile(basePath,[baseName,'.',eventsName,'.events.mat']);
end


if exist(filename,'file')
    eventstruct = load(filename);
    check_for_NWB = 0;
else
    check_for_NWB = 1;
end


if check_for_NWB
    if ~exist('eventsName','var')
        eventsName = [];        
    end
    eventstruct = loadEventsNWB(basePath, eventsName);
    if isempty(eventstruct)
       events = []; filename = [];
       return
    end
end


varsInFile = fieldnames(eventstruct);

if numel(varsInFile)==1
    events = eventstruct.(varsInFile{1});
else
    warning('Your .events.mat has multiple variables/structures in it... wtf.')
    events = eventstruct;
end

%Check that the events structure meets buzcode standards
[isEvents] = bz_isEvents(events);
switch isEvents
    case false
        warning('Your events structure does not meet buzcode standards. Sad.')
end



%% Function that loads NWB fields to an events structure
function eventstruct = loadEventsNWB(basePath, eventsName)
    % Check for NWB file since no .events.mat files were found
    d = dir(fullfile(basePath, '*.nwb'));
    if length(d) > 1 % If more than one .nwb files exist in the directory, select which one to load from
        warning('there is more than one .nwb file in this directory');
        [iNWBFile, ~] = listdlg('PromptString','Which NWB file would you like to load?',...
                                 'ListString',{d.name},'SelectionMode','single');
        d = d(iNWBFile);
    elseif length(d) == 0
        d = dir('*nwb');
        if isempty(d)
            warning('No .events.mat or NWB files present in this directory')
            eventstruct = [];
            return
        end
    end
    nwb_file = fullfile(d.folder, d.name); 
    
    if exist(nwb_file,'file') == 2
        nwb2 = nwbRead(nwb_file);
        
        % Check if the events field exists in the dataset
        events_key_exists = sum(ismember(keys(nwb2.processing),'events'));
        
        if events_key_exists
            all_event_keys = keys(nwb2.processing.get('events').nwbdatainterface);
        else
            disp('No events in this .nwb file')
            eventstruct = [];
            return
        end

        if isempty(all_event_keys)
            disp('No events in this .nwb file')
            eventstruct = [];
            return
        
        else
            % Check if a specific event was called to be loaded. If not, display a
            % pop-up list with the available events for selection
            if isempty(eventsName)
                [iEvent, ~] = listdlg('PromptString','Which event type would you like to load?',...
                                         'ListString',all_event_keys,'SelectionMode','single');
            else
                if sum(ismember(all_event_keys, eventsName))==1 % If the event exists, find its index
                    iEvent = find(ismember(all_event_keys, eventsName));
                else
                    disp('The selected event doesnt exist in this NWB file')
                    eventstruct = [];
                    return
                end
            end
            
            %% Get the selected event from the NWB file
            eventstruct = struct;
            
            event_module = nwb2.processing.get('events').nwbdatainterface.get(all_event_keys{iEvent}); % For convenience
            % Check if the timestamps are start-stop or just start
            % For information about the way this is coded, check:
            % https://nwb-schema.readthedocs.io/en/latest/format.html#sec-intervalseries
            first_Element = event_module.data.load(1);
            
            if ischar(first_Element) % This means I only have event start data
            	eventstruct.events.timestamps = event_module.timestamps.load;% neuroscope compatible matrix with 1 column - [starts] (in seconds)  
            else
                pointers = double(event_module.data.load)'; % positive values imply start, negative imply stop
                eventstruct.events.timestamps = [event_module.timestamps.load(pointers>0)' event_module.timestamps.load(pointers<0)'];% neuroscope compatible matrix with 2 columns - [starts stops] (in seconds)  
            end
            
            
            eventstruct.events.detectorinfo.detectorname       = 'N/A';% substructure with information about the detection method (fields below)
            eventstruct.events.detectorinfo.detectionparms     = 'N/A';% parameters used for detection  
            eventstruct.events.detectorinfo.detectiondate      = 'N/A';% date of detection
            eventstruct.events.detectorinfo.detectionintervals = 'N/A';% [start stop] pairs of intervals used for detection (optional)
%           eventstruct.events.detectorinfo.detectionchannel   = ;% channel used for detection (optional)  
        end
        
    else
        disp('No .events.mat or NWB files were found')
        eventstruct = [];
        return
    end
    
end
    

end