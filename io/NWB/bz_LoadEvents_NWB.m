function events = bz_LoadEvents_NWB(varargin)

%% Load events from NWB files following the Buzcode format

% Example calls:
% events = bz_LoadEvents_NWB;                                                    % Checks the current directory for a NWB file. Gives a list with the available events in the NWB file for the user to choose which one to load
% events = bz_LoadEvents_NWB('nwb_file', nwb_file);                              % Checks the specified NWB file. Gives a list with the available events in the NWB file for the user to choose which one to load
% events = bz_LoadEvents_NWB('nwb_file', nwb_file, 'PulseStim_5V_77777ms_LD12'); % Loads the specified event type from the NWB file
% events = bz_LoadEvents_NWB('PulseStim_5V_77777ms_LD12');                       % Checks the current directory for a NWB file. Loads the specified file


% The events are loaded straight from the NWB file. No extra files needed.

% This code is based on the bz_LoadEvents
% Konstantinos Nasiotis 2019


%%
    p = inputParser;
    addParameter(p,'nwb_file','',@isstr);
    
    if size(varargin,2)>2
        eventsName = varargin{3};
        varargin = {varargin{1,1:2}};
        parse(p,varargin{:})
        nwb_file = p.Results.nwb_file;
    elseif size(varargin,2)==1
        nwb_file   = [];
        eventsName = varargin{1};
    elseif size(varargin,2)==2
        parse(p,varargin{:})
        nwb_file = p.Results.nwb_file;
    elseif size(varargin,2)==0
        nwb_file   = [];
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
               error('could not find an nwb file..')
           end
       end
       nwb_file = fullfile(d.folder, d.name);
    end
    
    nwb2 = nwbRead(nwb_file);

    % Check if an events field exists in the dataset
    try
        events_exist_here = ~isempty(nwb2.stimulus_presentation);
        if ~events_exist_here
            disp('No events in this .nwb file')
            return
            
        else
            all_event_keys = keys(nwb2.stimulus_presentation);
            disp(' ')
            disp('The following event types are present in this dataset')
            disp('------------------------------------------------')
            for iEvent = 1:length(all_event_keys)
                disp(all_event_keys{iEvent})
            end
            disp(' ')
        end
    catch
        disp('No events in this .nwb file')
        return
    end
    
    
    % Check if a specific event was called to be loaded. If not, display a
    % pop-up list with the available events for selection
    if ~exist('eventsName','var')
        [iEvent, ~] = listdlg('PromptString','Which event type would you like to load?',...
                                 'ListString',all_event_keys,'SelectionMode','single');
    else
        if sum(ismember(all_event_keys, eventsName))>0
            iEvent = find(ismember(all_event_keys, eventsName));
        end
    end
    
    events = struct;

    events.timestamps                      = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).timestamps.load;% neuroscope compatible matrix with 1-2 columns - [starts stops] (in seconds)  
    events.detectorinfo.detectorname       = 'N/A';% substructure with information about the detection method (fields below)
    events.detectorinfo.detectionparms     = 'N/A';% parameters used for detection  
    events.detectorinfo.detectiondate      = 'N/A';% date of detection
    events.detectorinfo.detectionintervals = 'N/A';% [start stop] pairs of intervals used for detection (optional)
%         events.detectorinfo.detectionchannel   = ;% channel used for detection (optional)  

%         % Optional
%         events.amplitudes  = 1;% [Nx1 matrix]
%         events.frequencies = 1;% [Nx1 matrix]
%         events.durations   = 1;% [Nx1 matrix]


    % I add here the rest of the parameters that are saved on the nwb file
    events.starting_time_unit  = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).starting_time_unit;
    events.timestamps_interval = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).timestamps_interval;
    events.timestamps_unit     = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).timestamps_unit;
    events.comments            = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).comments;
    events.control             = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).control;
    events.control_description = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).control_description;
    events.data                = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).data;
    events.data_conversion     = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).data_conversion;
    events.data_resolution     = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).data_resolution;
    events.data_unit           = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).data_unit;
    events.description         = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).description;
    events.starting_time       = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).starting_time;
    events.starting_time_rate  = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).starting_time_rate;
    events.help                = nwb2.stimulus_presentation.get(all_event_keys{iEvent}).help;
    
    
end

