function [ EventDetectionReview ] = EventExplorer(basePath,events )
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
%OUTPUT
%   EventDetectionReview    results of DetectionReview - if called with an
%                           output, automatically runs detection review
%
%DLevenstein 2017
%% (For development)
% events = 'SlowWaves';
% basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%%
if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);

%% Select the Events
%Get the events structure, given the events input
if ~exist('events','var') %for no input - choose events
    [events,FO.eventsfilename] = bz_LoadEvents(basePath);
    if isempty(events)
        eventstype = 'none';
        eventsname = 'browse';
    else
        eventstype = 'events';
        [~,NAME,~] = fileparts(FO.eventsfilename);
        eventsname = NAME(length(baseName)+2:end-7); %baseName.eventsName.events
    end
elseif isstring(events) || ischar(events) %for eventsName - load the buzcode events.mad
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

%For events with start/stops, use the mean of the two as a marker. In the
%future replace with with start/stop markers
if size(exploreint,2)==2 %For events with start/stops, use the 
    exploreint = mean(exploreint,2);
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
    end
end

%Load EventExplorer metadata if it exists
FO.EEbuzcodefilename = fullfile(basePath,[baseName,'.EventExplorer.SessionMetadata.mat']);
if exist(FO.EEbuzcodefilename,'file')
    load(FO.EEbuzcodefilename)
    FO.FlagsAndComments.timepoint = EventExplorer.FlagsAndComments; 
    %This is an issue here - how to have general saved stuff (across
    %events) and event-specific saved stuff in events.mats.
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
    FO.lookatchannel = events.detectorinfo.detectionchannel;
catch
    FO.lookatchannel = inputdlg(['No events.detectorinfo.detectionchannel found in events.mat...',...
        'Which LFP channel would you like to look at?']);
    FO.lookatchannel = str2num(FO.lookatchannel{1});
end
    
%[ SleepState ] = bz_LoadStates(FO.basePath,'SleepState');
%FO.detectionints = SleepState.ints.NREMstate;
%FO.lookatchannel = events.detectorinfo.detectionparms.SWchannel;

FO.data.lfp = bz_GetLFP(FO.lookatchannel,'basepath',FO.basePath);
FO.data.spikes = bz_GetSpikes('basepath',FO.basePath);
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
FO.winsize = 8;
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
        'String',['(',num2str(length(FO.EventTimes)),' Total)'],...
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
            flaggedevents = FO.FlagsAndComments.(FO.viewmode).flags;
            %find the closest previous flagged event
            eventidx = interp1([0 flaggedevents],0:length(flaggedevents),...
                FO.currevent,'previous'); %set the current event to the closest previous flagged event
            FO.currevent = flaggedevents(eventidx+1); %next flagged event
        otherwise
            switch FO.viewmode
                case 'timepoint'
                    FO.currevent=FO.currevent+0.75.*FO.winsize;
                otherwise
                    FO.currevent=FO.currevent+1;
            end
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
            FO.currevent=randi(length(FO.EventTimes));
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

    switch FO.viewmode
        case 'timepoint'
            %THis is going to be buggy with multiple timepoint.flags for
            %timepoint.comments......
        otherwise 
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
            
            try %also flag the timepoint for events - going to be issue with comments
                FO.FlagsAndComments.timepoint.flags(end+1) = FO.EventTimes(FO.currevent);
                FO.FlagsAndComments.timepoint.comments{end+1} = [];
            catch
                FO.FlagsAndComments.timepoint.flags = FO.EventTimes(FO.currevent);
                FO.FlagsAndComments.timepoint.comments{1} = [];
            end
    end
    [FO.FlagsAndComments.timepoint.flags,I] = sort(FO.FlagsAndComments.timepoint.flags);
    FO.FlagsAndComments.timepoint.comments=FO.FlagsAndComments.timepoint.comments(I);
    FO.FlagsAndComments.(FO.viewmode).flags = sort(FO.FlagsAndComments.(FO.viewmode).flags);
    guidata(FO.fig, FO);
end

function AddUserComment(obj,event)
    FO = guidata(obj);
    usercomment = get(obj,'String');
    switch FO.viewmode
        case 'timepoint' %If browsing timepoints, need to record time and comment
            try
                FO.FlagsAndComments.timepoint.comments{end+1} = usercomment;
                FO.FlagsAndComments.timepoint.flags(end+1) = FO.currevent;
            catch
                FO.FlagsAndComments.timepoint.comments{1} = usercomment;
                FO.FlagsAndComments.timepoint.flags = FO.currevent;
            end
            set(FO.flageventbutton,'String','Unflag')
        otherwise %if events/misses/FA, can just use current event number
            FO.FlagsAndComments.(FO.viewmode).comments{FO.currevent} = usercomment;
            try %also comment the time point (this is redundant)
                FO.FlagsAndComments.timepoint.comments{end+1} = usercomment;
                FO.FlagsAndComments.timepoint.flags(end+1) = FO.EventTimes(FO.currevent);
            catch
                FO.FlagsAndComments.timepoint.comments{1} = usercomment;
                FO.FlagsAndComments.timepoint.flags = FO.EventTimes(FO.currevent);
            end
    end
    [FO.FlagsAndComments.timepoint.flags,I] = sort(FO.FlagsAndComments.timepoint.flags);
    FO.FlagsAndComments.timepoint.comments=FO.FlagsAndComments.timepoint.comments(I);
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

function ChangeLFPChan(obj,event)
    FO = guidata(obj);
    FO.lookatchannel=str2num(obj.String);
    FO.data.lfp = bz_GetLFP(FO.lookatchannel,'basepath',FO.basePath);
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
                    eventsfile.(FO.EventName).EventExplorer.FlagsAndComments = FO.FlagsAndComments;
                    save(FO.eventsfilename,'-struct','eventsfile',FO.EventName,'-append')
                catch
                    warndlg({' Save failed... ',[FO.eventsfilename,' may not ',...
                        'contain a structure titled ',FO.EventName,'.'],...
                        'Or you may not have sudo priviliges...?'},'Oh No!')
                end
        end
        %Save the General EventExplorer metadata file
        EventExplorer.FlagsAndComments = FO.FlagsAndComments.timepoint;
        save(FO.EEbuzcodefilename,'EventExplorer')
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