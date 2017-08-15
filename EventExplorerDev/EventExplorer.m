function [ EEoutput ] = EventExplorer(basePath,events )
%EventExplorer is a GUI tool for exploring buzcode events and states files.
%
%INPUT
%   events      Name of events in a string (i.e. 'SlowWaves') or buzcode 
%               events structure
%   basePath    default: pwd
%
%DLevenstein 2017
%%
%events = 'SlowWaves';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%%
if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);


%% Select the Events
%eventsfile = fullfile(basePath,[baseName,'.',eventsname,'.events.mat']);
if ~exist('events','var')
    [events,FO.eventsfilename] = bz_LoadEvents(basePath);
    eventstype = 'events';
    [~,NAME,~] = fileparts(FO.eventsfilename);
    eventsname = NAME(length(baseName)+2:end-7); %baseName.eventsName.events
elseif isstring(events) || ischar(events)
    eventsname = events;
    [events,FO.eventsfilename] = bz_LoadEvents(basePath,events);
    eventstype = 'events';
    
    if isempty(events)
        display('no events.mat by that name, trying to find a states.mat')
        [events,FO.eventsfilename] = bz_LoadStates(basePath,events);
        eventstype = 'states';
    end
else
    eventstype = 'events'; %Need a way to check type of a struct?
    eventsname = inputname(2);
end

switch eventstype
    case 'events'
        exploreint = events.timestamps;
        exploreintname = events.detectorinfo.detectorname;

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

end
%% 
FO.baseName = baseName;
FO.EventTimes = exploreint;
FO.EventName = eventsname;
FO.basePath = basePath;

%Load EventExplorer data from events file 
REVIEWDONE = false;
if isfield(events,'EventExplorer')
    if isfield(events.EventExplorer,'FlagsAndComments')
        FO.FlagsAndComments = events.EventExplorer.FlagsAndComments;
    end
    if isfield(events.EventExplorer,'DetectionReview')
        REVIEWDONE=true;
        FO.DetectionReview = events.EventExplorer.DetectionReview;
    elseif isfield(events,'EventReview') %For legacy SWDetection - remove in next few days (8/15)
        FO.DetectionReview = events.EventReview;
    end
end
%% Load The Data, eh?
%could put to function: EE_Initiate
try
    FO.detectionints = events.detectorinfo.detectionintervals;
catch
    display('No events.detectorinfo.detectionintervals found in events.mat... using all time')
    FO.detectionints = [0 Inf];
end
try
    FO.detectionchannel = events.detectorinfo.detectionchannel;
catch
    try %For legacy SWDetection - remove in next few days (8/15)
        FO.detectionchannel = events.detectorinfo.detectionparms.SWchannel;
    catch
    FO.detectionchannel = inputdlg(['No events.detectorinfo.detectionchannel found in events.mat...',...
        'Which LFP channel would you like to look at?']);
    FO.detectionchannel = str2num(FO.detectionchannel{1});
    end
end
    
%[ SleepState ] = bz_LoadStates(FO.basePath,'SleepState');
%FO.detectionints = SleepState.ints.NREMstate;
%FO.detectionchannel = events.detectorinfo.detectionparms.SWchannel;

FO.data.lfp = bz_GetLFP(FO.detectionchannel,'basepath',FO.basePath);
FO.data.spikes = bz_GetSpikes('basepath',FO.basePath);
%% Set up the EventExplorer Window
%Position for the main interface
posvar = get(0,'Screensize');
posvar(1) = 20;
posvar(2) = 20;
posvar(3) = posvar(3)-100;
posvar(4) = posvar(4)-100;

oldfig = findobj('tag','EventExplorerMaster'); close(oldfig);
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
FO.winsize = 8;
FO.currevent = 1;
FO.viewmode = 'events';


%Text of hotkey definitions for user guidance

%Set up the navigation panel
FO.NavPanel = uipanel('FontSize',12,...
        'Position',[.65 .22 0.25 0.15]);    
    thisevent  = uicontrol('Parent',FO.NavPanel,'Style','edit',...
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
    r1 = uicontrol(FO.eventtypeselection,'Style',...
                      'radiobutton',...
                      'String','events',...
                      'Position',[10 70 75 30]);
    if REVIEWDONE
    r2 = uicontrol(FO.eventtypeselection,'Style','radiobutton',...
                      'String','misses',...
                      'Position',[10 40 75 30]);
    r3 = uicontrol(FO.eventtypeselection,'Style','radiobutton',...
                      'String','FAs',...
                      'Position',[10 10 75 30]);
    end
    randbtn = uicontrol('Parent',FO.eventtypeselection,...
        'Position',[130 40 150 40],'String','Run Detection Review',...
         'Callback',@DetectionReview);

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
        
%Store the data in the figure - do this at the end of each function?
guidata(FO.fig, FO);
EventVewPlot;
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
            flaggedevents = FO.FlagsAndComments.(FO.viewmode).flags;
            %find the closest previous flagged event
            eventidx = interp1([0 flaggedevents],0:length(flaggedevents),...
                FO.currevent,'previous'); %set the current event to the closest previous flagged event
            FO.currevent = flaggedevents(eventidx+1); %next flagged event
        otherwise
            FO.currevent=FO.currevent+1;
    end
    guidata(FO.fig, FO);
    EventVewPlot;
end

function PrevEvent(obj,eventdata)
    FO = guidata(obj); 
    switch FO.showflagged
        case true
            flaggedevents = FO.FlagsAndComments.(FO.viewmode).flags;
            %find the closest previous flagged event
            eventidx = interp1([flaggedevents max([flaggedevents(end) FO.currevent])+1],...
                1:length(flaggedevents)+1,...
                FO.currevent,'next'); %set the current event to the closest next flagged event
            FO.currevent = flaggedevents(eventidx-1); %previous flagged event
        otherwise
            FO.currevent=max(FO.currevent-1,1);
    end
    guidata(FO.fig, FO);
    EventVewPlot;
end

function RandEvent(obj,eventdata)
    FO = guidata(obj);
    switch FO.viewmode
        case 'events'
            FO.currevent=randi(length(FO.EventTimes));
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
    try  %This is to deal with case where FO.FlagsAndComments.(FO.viewmode).flags hasn't been made yet... do better
        [isflagged,flagidx] = ismember(FO.currevent,FO.FlagsAndComments.(FO.viewmode).flags);
        switch isflagged
            case true
                FO.currevent,FO.FlagsAndComments.(FO.viewmode).flags(flagidx)=[];
                set(FO.flageventbutton,'String','Flag')
            case false
                FO.FlagsAndComments.(FO.viewmode).flags(end+1) = FO.currevent;
                set(FO.flageventbutton,'String','Unflag')
        end
    catch %If flags haven't yet been added to FO for this viewmode
        FO.FlagsAndComments.(FO.viewmode).flags = FO.currevent;
        set(FO.flageventbutton,'String','Unflag')
    end
    FO.FlagsAndComments.(FO.viewmode).flags = sort(FO.FlagsAndComments.(FO.viewmode).flags);
    guidata(FO.fig, FO);
end

function AddUserComment(obj,event)
    FO = guidata(obj);
    usercomment = get(obj,'String');
    FO.FlagsAndComments.(FO.viewmode).comments{FO.currevent} = usercomment;
    guidata(FO.fig, FO);
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
                    eventsfile.(FO.EventName).EventExplorer.FlagsAndComments = FO.FlagsAndComments;
                    save(FO.eventsfilename,'-struct','eventsfile',FO.EventName,'-append')
                catch
                    warndlg({' Save failed... ',[FO.eventsfilename,' may not ',...
                        'contain a structure titled ',FO.EventName,'.'],...
                        'Or you may not have sudo priviliges...?'},'Oh No!')
                end
        end
    end
delete(FO.fig)
end
