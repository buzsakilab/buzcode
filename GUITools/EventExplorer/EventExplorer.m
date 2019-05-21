function [ EventDetectionReview ] = EventExplorer(basePath,events,varargin )
%EventExplorer is a GUI tool for exploring buzcode events and states files.
%
%INPUT
%   events      eventsName string (i.e. 'SlowWaves'), 
%               to load file: basePath/baseName.eventsName.events.mat
%                   -or- 
%               buzcode events structure (see also buzcode wiki)
%               Required fields: 
%                   eventsName.timestamps  [numevents x 1] or [numevents x 2]
%               Recommended fields:
%                   eventsName.detectorinfo.detectionintervals
%                   eventsName.detectorinfo.detectionchannel
%                   
%   basePath    basePath that holds baseName.lfp (default: pwd)
%
%(options)
%   'useSpikes' true,false,'MUA'
%   'noPrompts' default: true
%
%OUTPUT
%   EventDetectionReview    results of DetectionReview - if called with an
%                           output, automatically runs detection review
%
%DLevenstein 2017
%% (For development)
% events = 'SlowWaves';
% basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%%
p = inputParser;
addParameter(p,'useSpikes',false);
addParameter(p,'noPrompts',true,@islogical);

parse(p,varargin{:})
noPrompts = p.Results.noPrompts;
useSpikes = p.Results.useSpikes;



%%
if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);

%% Select the Events (this should go into an internal function: LoadEvents
%Get the events structure, given the events input
if ~exist('events','var') %for no input - choose events
    [events,FO.eventsfilename] = bz_LoadEvents(basePath);
    if isempty(events) %Needed: way to get to eventstype - none...
        eventstype = 'none';
        eventsname = 'browse';
    else
        eventstype = 'events';
        [~,NAME,~] = fileparts(FO.eventsfilename);
        eventsname = NAME(length(baseName)+2:end-7); %baseName.eventsName.events
    end
elseif isstring(events) || ischar(events) %for eventsName - load the buzcode events.mat
    eventsname = events;
    [events,FO.eventsfilename] = bz_LoadEvents(basePath,events);
    eventstype = 'events';
    if isempty(events)
        display('no events.mat by that name, trying to find a states.mat')
        [events,FO.eventsfilename] = bz_LoadStates(basePath,events);
        eventstype = 'states';
    end
else %for a buzcode structure

    eventstype = 'events'; %Need a way to check type of a struct?
    eventsname = inputname(2);
end

%Get the right info out of the events structure, given its type
%Sets the view mode to either 'events' (viewing pre-marked events) or
%'timepoint' = browsing time points that have been previously
%flagged/commented
switch eventstype
    case 'events'
        exploreint = events.timestamps;
        %exploreintname = events.detectorinfo.detectorname;
        FO.viewmode = 'events';
    case 'states'
        intnames = [fieldnames(events.ints)]; %remove second term here
        [s,v] = listdlg('PromptString','Which interval type to explore?',...
            'SelectionMode','single','ListString',intnames);
        exploreintname = intnames{s};
        try
        exploreint = events.ints.(exploreintname);
        catch
            exploreint = events.(exploreintname);
        end
        FO.viewmode = 'events';
    case 'none'
        exploreint = [];
        FO.viewmode = 'timepoint';
end

%% 
FO.baseName = baseName;
FO.EventName = eventsname;
FO.basePath = basePath;

%NOTE: F0.EventTimes should either be Nevents x 1 or Nevents x 2 (starts/stops).
%If (starts/stops), flags/etc will be associated with the start time
FO.EventTimes = exploreint;

%Load any EventExplorer data from events file (events FaC)
REVIEWDONE = false;
if isfield(events,'EventExplorer')
    if isfield(events.EventExplorer,'FlagsAndComments')
        FO.FlagsAndComments.events = events.EventExplorer.FlagsAndComments;
    else
        FO.FlagsAndComments.events.flags = [];
        FO.FlagsAndComments.events.comments = {};
    end
    if isfield(events.EventExplorer,'DetectionReview')
        REVIEWDONE=true;
        FO.DetectionReview = events.EventExplorer.DetectionReview;
    end
else
	FO.FlagsAndComments.events.flags = [];
	FO.FlagsAndComments.events.comments = {};
end

%Load EventExplorer metadata if it exists (timestamps FaC)
FO.EEbuzcodefilename = fullfile(basePath,[baseName,'.EventExplorer.SessionMetadata.mat']);
if exist(FO.EEbuzcodefilename,'file')
    load(FO.EEbuzcodefilename)
    FO.FlagsAndComments.timepoint = EventExplorer.FlagsAndComments; 
else
	FO.FlagsAndComments.timepoint.flags = [];
	FO.FlagsAndComments.timepoint.comments = {};
end

%Align any preexisting FaC (timestamps and events)
FO.FlagsAndComments = MergeFlagsComments(FO.FlagsAndComments,FO.EventTimes);
%% Load The Data, eh?
%could put to function: EE_Initiate

%Get the intervals,channel used for detection out of the events.mat file
try
    FO.detectionints = events.detectorinfo.detectionintervals;
catch
    display('No events.detectorinfo.detectionintervals found in events.mat... using all time')
    FO.detectionints = [0 Inf];
end
try
    FO.lookatchannel = events.detectorinfo.detectionchannel;
catch
    FO.lookatchannel = inputdlg(['No events.detectorinfo.detectionchannel found in events.mat...',...
        'Which LFP channel would you like to look at?']);
    if isempty(FO.lookatchannel)
        return
    end
    FO.lookatchannel = str2num(FO.lookatchannel{1});
end
   
%Load the LFP and spikes
FO.data.lfp = bz_GetLFP(FO.lookatchannel,'basepath',FO.basePath,'noPrompts',true);
switch useSpikes
    case true
        FO.data.spikes = bz_GetSpikes('basepath',FO.basePath);
    case false
        FO.data.spikes = [];
    case 'MUA'
        MUA = MUAfromDat(basePath,'saveMat',true);
        FO.data.spikes = MUA.peaks;
end
%% Set up the EventExplorer Window
%Position for the main interface
posvar = get(0,'Screensize');
posvar(1) = 20;
posvar(2) = 20;
posvar(3) = posvar(3)-100;
posvar(4) = posvar(4)-100;

oldfig = findobj('tag','EventExplorerMaster'); 
try close(oldfig); catch;  delete(oldfig); end;
%Start the figure
FO.fig = figure('KeyPressFcn', {@KeyDefinitions},'Position', posvar);
set(FO.fig, 'numbertitle', 'off', 'name', ['Recording: ', FO.baseName,'. Events: ',FO.EventName]);
set(FO.fig, 'Tag', 'EventExplorerMaster','menubar', 'none');
set(FO.fig, 'CloseRequestFcn', {@CloseDialog});
set(FO.fig,'WindowButtonDownFcn', {@MouseClick});

%From StateEditor - anything else here needed?
% set(FO.fig,'WindowButtonDownFcn', {@MouseClick}, 'WindowButtonUpFcn', {@unMouseClick}, 'Units', 'normalized');
% set(FO.fig, 'WindowButtonMotionFcn', {@Nothing}, 'WindowScrollWheelFcn', {@MouseScroll});

FO.currentuseraction = 'none';

%Set up the view window
FO.viewwin = subplot(3,1,2,'ButtonDownFcn',@MouseClick);
%Event Selection panel

FO.scaleLFP = 1;
FO.winsize = 8; %8second window
FO.currevent = 1;


%Text of hotkey definitions for user guidance


%The SignalBox
FO.SignalPanel = uipanel('FontSize',12,...
        'Position',[.65 .7 0.25 0.15]); 
    editLFPchan  = uicontrol('Parent',FO.SignalPanel,'Style','edit',...
        'Position',[230 70 50 25],'String',num2str(FO.lookatchannel),...
        'Callback',@ChangeLFPChan);
    winsizetext = uicontrol('Parent',FO.SignalPanel,...
        'Position',[150 70 80 18],'style','text',...
        'string','LFP Channel:','HorizontalAlignment','left'); 
    loadbehavior = uicontrol('Parent',FO.SignalPanel,...
        'Position',[230 90 50 25],'String','Load Behavior',...
         'Callback',@LoadBeh);


%Set up the navigation panel
FO.NavPanel = uipanel('FontSize',12,...
        'Position',[.65 .22 0.25 0.15]);    
    FO.thiseventdisplay  = uicontrol('Parent',FO.NavPanel,'Style','edit',...
        'Position',[60 70 50 25],'String',num2str(FO.currevent),...
        'Callback',@GoToEvent);
    thiseventtext = uicontrol('Parent',FO.NavPanel,...
        'Position',[20 70 40 18],'style','text',...
        'string','Event #','HorizontalAlignment','left'); 
    nextbtn = uicontrol('Parent',FO.NavPanel,...
        'Position',[110 20 60 40],'String','->',...
         'Callback',@NextEvent);
    prevbtn = uicontrol('Parent',FO.NavPanel,...
        'Position',[40 20 60 40],'String','<-',...
         'Callback',@PrevEvent);
    randbtn = uicontrol('Parent',FO.NavPanel,...
        'Position',[200 20 100 40],'String','Random',...
         'Callback',@RandEvent);
    editwinsize  = uicontrol('Parent',FO.NavPanel,'Style','edit',...
        'Position',[230 70 50 25],'String',num2str(FO.winsize),...
        'Callback',@EditWinSize);
    winsizetext = uicontrol('Parent',FO.NavPanel,...
        'Position',[150 70 80 18],'style','text',...
        'string','Window Size (s):','HorizontalAlignment','left'); 

%The Reviewed Event Selection Panel
FO.eventtypeselection = uibuttongroup('Position',[0.65,0.05,0.25,0.15],'Visible','on',...
    'SelectionChangedFcn',@(bg,event) EventTypeSelector(bg,event));
    eventselectiontext = uicontrol(FO.eventtypeselection,'Style','radiobutton',...
                      'String','events','Position',[10 70 75 30]);
    FO.missbutton = uicontrol(FO.eventtypeselection,'Style','radiobutton',...
                      'String','misses','Visible','off',...
                      'Position',[10 40 75 30]);
    FO.FAbutton = uicontrol(FO.eventtypeselection,'Style','radiobutton',...
                      'String','FAs','Visible','off',...
                      'Position',[10 10 75 30]);
    FO.eventcounttxt = uicontrol(FO.eventtypeselection,'Style','text',...
        'String',['(',num2str(size(FO.EventTimes,1)),' Total)'],...
        'Position',[75 70 80 22]);
    
    FO.missperctxt = uicontrol(FO.eventtypeselection,'Style','text',...
        'String','',...
        'Position',[75 40 80 22]);
    FO.FAperctxt = uicontrol(FO.eventtypeselection,'Style','text',...
        'String','',...
        'Position',[75 10 80 22]);
    if REVIEWDONE
        set(FO.missbutton,'Visible','on');
        set(FO.missperctxt,'String',['(Est ',num2str(round(FO.DetectionReview.estMissperc.*100,2)),'%)']);
        set(FO.FAbutton,'Visible','on');
        set(FO.FAperctxt,'String',['(Est ',num2str(round(FO.DetectionReview.estFAperc.*100,2)),'%)']);
    end
    rundetectionbtn = uicontrol('Parent',FO.eventtypeselection,...
        'Position',[165 40 150 40],'String','Run Detection Review',...
         'Callback',@RunDetectionReview);

%The Comment/Flag Panel
FO.showflagged = 'false';  %state of the Flagged Only checker
FO.CommentFlagPanel = uipanel('FontSize',12,...
        'Position',[.1 .05 0.5 0.3]);
    FO.flageventbutton = uicontrol('Parent',FO.CommentFlagPanel,...
        'Position',[550 20 100 40],'String','Flag',...
         'Callback',@FlagEvent);
    FO.eventcomment = uicontrol('Parent',FO.CommentFlagPanel,...
        'Style','edit','Max',5,'String','Event Comments',...
        'HorizontalAlignment','left','Position',[20 60 630 160],...
        'Callback',@AddUserComment);
    flaggedonly = uicontrol('Parent',FO.CommentFlagPanel,'Style','checkbox',...
        'Position',[430 20 20 40],...
        'Callback',@ShowFlagged);
    showflagtext = uicontrol('Parent',FO.CommentFlagPanel,...
        'Position',[450 20 100 30],'style','text',...
        'string','Browse Flagged Events Only','HorizontalAlignment','left'); 
        
    
%Store the data in the figure guidata
guidata(FO.fig, FO);
EventVewPlot;

%If there's an output - run DetectionReview
if nargout >0
    DetectionReview(FO.fig);
    obj = findobj('tag','EventExplorerMaster');  FO = guidata(obj); %get the results from the guidata
    EventDetectionReview = FO.DetectionReview;
end

end %Gen function end. Below are callback definitions

%% %%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALLBACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function KeyDefinitions(f, e)
    obj = findobj('tag','EventExplorerMaster');  FO = guidata(obj); 
    switch e.Key
        case 'uparrow'
            FO.scaleLFP = FO.scaleLFP+0.1;
            guidata(FO.fig, FO); 
            EventVewPlot;
        case 'downarrow'
            FO.scaleLFP = FO.scaleLFP-0.1;
            guidata(FO.fig, FO); 
            EventVewPlot;
        case 'rightarrow';  NextEvent(obj);
        case 'leftarrow';   PrevEvent(obj)
    end
end

function NextEvent(obj,eventdata)
    FO = guidata(obj); 
    switch FO.showflagged
        case true
            %For Browsing flagged events only
            flaggedevents = FO.FlagsAndComments.(FO.viewmode).flags;
            if ~isempty(flaggedevents)
                %find the closest previous flagged event
                eventidx = interp1([0 flaggedevents],0:length(flaggedevents),...
                    FO.currevent,'previous'); %set the current event to the closest previous flagged event
                FO.currevent = flaggedevents(min(eventidx+1,length(flaggedevents))); %next flagged event
            end
        otherwise
            switch FO.viewmode
                case 'timepoint'
                    %For timepoint mode: jump 75% of windowsize
                    FO.currevent=FO.currevent+0.75.*FO.winsize;
                otherwise
                    FO.currevent=min(FO.currevent+1,size(FO.EventTimes,1));
            end
    end
    guidata(FO.fig, FO);
    EventVewPlot;
end

function PrevEvent(obj,eventdata)
    FO = guidata(obj); 
    switch FO.showflagged
        case true
            %For Browsing flagged events only
            flaggedevents = FO.FlagsAndComments.(FO.viewmode).flags;
            if ~isempty(flaggedevents)
                %find the closest previous flagged event
                eventidx = interp1([flaggedevents max([flaggedevents(end) FO.currevent])+1],...
                    1:length(flaggedevents)+1,...
                    FO.currevent,'next'); %set the current event to the closest next flagged event
                FO.currevent = flaggedevents(max(eventidx-1,1)); %previous flagged event
            end
        otherwise
            switch FO.viewmode
                case 'timepoint'
                    FO.currevent=FO.currevent-0.75.*FO.winsize;
                otherwise
                    FO.currevent=max(FO.currevent-1,1);
            end
    end
    guidata(FO.fig, FO);
    EventVewPlot;
end

function RandEvent(obj,eventdata)
    FO = guidata(obj);
    switch FO.viewmode
        case 'events'
            FO.currevent=randi(size(FO.EventTimes,1));
        case 'timepoint'
            FO.currevent=randsample(FO.data.lfp.timestamps);
        otherwise
            display(['Random does not yet work for viewmode "',FO.viewmode,'"'])
    end
    guidata(FO.fig, FO);
    EventVewPlot;
end

function GoToEvent(obj,eventdata)
    FO = guidata(obj);
    switch FO.viewmode
        case 'events'
            FO.currevent=str2num(obj.String);
        case 'timepoint'
            FO.currevent=str2num(obj.String);
        otherwise
            display(['GoTo does not yet work for viewmode "',FO.viewmode,'"'])
    end
    guidata(FO.fig, FO);
    EventVewPlot;
end

function EditWinSize(obj,eventdata)
    FO = guidata(obj);
    FO.winsize=str2num(obj.String);
    guidata(FO.fig, FO);
    EventVewPlot;
end

function EventTypeSelector(obj,event)
    FO = guidata(obj);
    FO.viewmode = event.NewValue.String;
    FO.currevent = 1;
    guidata(FO.fig, FO);
    EventVewPlot;
end

function FlagEvent(obj,event)
    FO = guidata(obj); 
    usercomment = FO.eventcomment.String;
    if strcmp(usercomment,'Event Comments')
        usercomment = [];
    end
    
    %Has it already been flagged?
    [isflagged,flagidx] = ismember(FO.currevent,FO.FlagsAndComments.(FO.viewmode).flags);
    switch isflagged
        case false
            %Flag it! (and the time point)
            FO.FlagsAndComments.(FO.viewmode).flags(end+1) = FO.currevent;
            FO.FlagsAndComments.(FO.viewmode).comments{end+1} = usercomment;
            set(FO.flageventbutton,'String','Unflag')
        case true
            %Remove the Flag/Comment
            FO.FlagsAndComments.(FO.viewmode).flags(flagidx)=[];
            FO.FlagsAndComments.(FO.viewmode).comments(flagidx) = [];
            set(FO.flageventbutton,'String','Flag')
    end
    
    FO.FlagsAndComments = MergeFlagsComments(FO.FlagsAndComments,FO.EventTimes);
    guidata(FO.fig, FO);
end

function AddUserComment(obj,event)
    FO = guidata(obj);
    usercomment = get(obj,'String');
    
    %If it's flagged already, update the comment. If not, do nothing
    [isflagged,flagidx] = ismember(FO.currevent,FO.FlagsAndComments.(FO.viewmode).flags);
    if isflagged
        FO.FlagsAndComments.(FO.viewmode).comments{flagidx} = usercomment;
        
        %Find the timepoint
        if ~strcmp(FO.viewmode,'timepoint')
            [~,flagtimeidx] = ismember(FO.EventTimes(FO.currevent,1),FO.FlagsAndComments.timepoint.flags);
            FO.FlagsAndComments.timepoint.comments{flagtimeidx} = usercomment;
        end

    end
    
    guidata(FO.fig, FO);
end

function FlagsAndComments = MergeFlagsComments(FlagsAndComments,EventTimes)
    %Merges the timepoint/event comments/flags to maintain continuity.
    %Any flagged timepoints that are very near an event are moved to that
    %event time, and the event is flagged/commented. The timepoints of
    %any flagged events are flagged.
    %NOTE: this could be done more efficiently with uniquetol or
    %ismembertol
    timetol = 0.02; %within 20ms
    
    %Make sure FlagsAndComments has both events and timepoints
    
    
    %Bring the flagged timepoints to the flagged events
    flaggedtimepoints = FlagsAndComments.timepoint.flags;
    flaggedeventtimes = EventTimes(FlagsAndComments.events.flags,1);
    for tp = 1:length(flaggedtimepoints)
        durtoflag = abs(flaggedtimepoints(tp)-EventTimes(:,1));

        %If there's an event within tolerance of the time point 
        if any(durtoflag<timetol)
            %move the time point to the close event    
            FlagsAndComments.timepoint.flags(tp) = EventTimes(durtoflag<timetol,1);
            %Is its not already in the list of flagged events?
            if ~ismember(FlagsAndComments.timepoint.flags(tp),flaggedeventtimes)
                %Which event is it? Add it!
                FlagsAndComments.events.flags(end+1)=...
                    find(ismember(EventTimes(:,1),FlagsAndComments.timepoint.flags(tp)));
                FlagsAndComments.events.comments{end+1}=...
                    FlagsAndComments.timepoint.comments{tp};
            end
                
        end
    end
    
    %Bring the flagged events to the flagged timepoints
    flaggedeventtimes = EventTimes(FlagsAndComments.events.flags,1);
    for et = 1:length(flaggedeventtimes)
        durtoflag = abs(flaggedeventtimes(et)-FlagsAndComments.timepoint.flags);
        
        if any(durtoflag<timetol)
            %If there's a timpoint within tolerance, move the timepoint to the event time
            FlagsAndComments.timepoint.flags(durtoflag<timetol)=flaggedeventtimes(et);
        else
            %If not, add the timepoint of the event.
            FlagsAndComments.timepoint.flags(end+1)=flaggedeventtimes(et);
            FlagsAndComments.timepoint.comments{end+1}=FlagsAndComments.events.comments{et};
        end
    end
           
    %Sort everything and return it
    [FlagsAndComments.timepoint.flags,I] = sort(FlagsAndComments.timepoint.flags);
    FlagsAndComments.timepoint.comments=FlagsAndComments.timepoint.comments(I);
    
    [FlagsAndComments.events.flags,I] = sort(FlagsAndComments.events.flags);
    FlagsAndComments.events.comments=FlagsAndComments.events.comments(I);
  
    %Merge any duplicates...
%     if any(diff(FlagsAndComments.events.flags)==0)
%         FlagsAndComments.events.flags(diff(FlagsAndComments.events.flags)==0)
%         any(diff(FlagsAndComments.events.flags)==0)
%     end
end

function ShowFlagged(obj,event) 
    FO = guidata(obj);
    if (get(obj,'Value') == get(obj,'Max'))
        FO.showflagged = true;
    else
        FO.showflagged = false;
    end
    guidata(FO.fig, FO);
end

function ChangeLFPChan(obj,event)
    FO = guidata(obj);
    FO.lookatchannel=str2num(obj.String);
    FO.data.lfp = bz_GetLFP(FO.lookatchannel,'basepath',FO.basePath);
    guidata(FO.fig, FO);
    EventVewPlot;
end

function LoadBeh(obj,event)
    FO = guidata(obj);
    FO.lookatchannel=str2num(obj.String);
    if isfield(FO,'behavior')
        numbeh = length(FO.behavior);
        FO.behavior{numbeh+1} = bz_LoadBehavior(FO.basePath);
    else
        numbeh = 0;
        FO.behavior{1} = bz_LoadBehavior(FO.basePath);
    end
    if ~isfield(FO.behavior{numbeh+1},'data')
        hasfields = fieldnames(FO.behavior{numbeh+1});
        [s,v] = listdlg('PromptString','Which field is your data?',...
                     'ListString',hasfields,'SelectionMode','single');
        FO.behavior{numbeh+1}.data = FO.behavior{numbeh+1}.(hasfields{s});
    end
    guidata(FO.fig, FO);
    EventVewPlot;
end

function CloseDialog(obj,event)
FO = guidata(obj);
    if isfield(FO,'eventsfilename') && isfield(FO,'FlagsAndComments')
        %Need to check here if no changes were made using ISEQUAL(A,B) no
        %prompt if the saved stuff is same as FO.
        button = questdlg(['Would you like to save flags/comments to ',...
            FO.eventsfilename,'?'],'Good Bye.');
        switch button
            case 'Yes'
                %Load the events file, add the field, save the events file
                try %Only do this if the correct named structure lives in the file
                    eventsfile = load(FO.eventsfilename,FO.EventName);
                    eventsfile.(FO.EventName).EventExplorer.FlagsAndComments = FO.FlagsAndComments.events;
                    save(FO.eventsfilename,'-struct','eventsfile',FO.EventName,'-append')
                catch
                    warndlg({' Save failed... ',[FO.eventsfilename,' may not ',...
                        'contain a structure titled ',FO.EventName,'.'],...
                        'Or you may not have sudo priviliges...?'},'Oh No!')
                end
        end
        
        %Save the General EventExplorer metadata file
        if isfield(FO,'FlagsAndComments')
            EventExplorer.FlagsAndComments = FO.FlagsAndComments.timepoint;
            try
                save(FO.EEbuzcodefilename,'EventExplorer')
            catch
                disp('unable to save Flags and Comments')
            end
        end
    end
delete(FO.fig)
end

function RunDetectionReview(obj,event)
    DetectionReview(obj); %Run
    %Get the results and Update miss/etc tickers
    obj = findobj('tag','EventExplorerMaster');  FO = guidata(obj);
        set(FO.missbutton,'Visible','on');
        set(FO.missperctxt,'String',['(Est ',num2str(round(FO.DetectionReview.estMissperc,2)),'%)']);
        set(FO.FAbutton,'Visible','on');
        set(FO.FAperctxt,'String',['(Est ',num2str(round(FO.DetectionReview.estFAperc,2)),'%)']);
end