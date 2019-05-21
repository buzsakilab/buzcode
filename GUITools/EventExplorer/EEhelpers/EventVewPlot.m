function [ viewinfo ] = EventVewPlot
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
obj = findobj('tag','EventExplorerMaster');  FO = guidata(obj); 

%Two options of current timepoint: Event Number, arbitraty timepoint
%FA/misses for DetectionReview
switch FO.viewmode
    case 'timepoint' 
        timepoint = FO.currevent;
    case 'events'
        %Always use the start. Could update to midpoint if needed
        timepoint = FO.EventTimes(FO.currevent,1); 
    case 'FAs'
        timepoint = FO.DetectionReview.falsealarm(FO.currevent);
    case 'misses'
        timepoint = FO.DetectionReview.miss(FO.currevent);
end

%Get the view window and events inside it
thiseventwin = timepoint+FO.winsize.*[-0.5 0.5];

%Plot
%viewwin = subplot(3,1,2,'ButtonDownFcn',@MouseClick);
%set(gca,'ButtonDownFcn', @MouseClick)
hold(FO.viewwin,'off')
%Plot The LFP
ywinrange = bz_MultiLFPPlot(FO.data.lfp,'timewin',thiseventwin,'spikes',FO.data.spikes,...
    'axhandle',FO.viewwin,'scaleLFP',FO.scaleLFP)
hold on

%Plot the events
if size(FO.EventTimes,2) == 1 %Single timepoint events
    inwinevents = FO.EventTimes>=thiseventwin(1) & FO.EventTimes<=thiseventwin(2);
    inwineventtimes = FO.EventTimes(inwinevents);
    plot(FO.viewwin,inwineventtimes,zeros(size(inwineventtimes)),...
        'o','color',[0 0.6 0],'markersize',5,'linewidth',3)

elseif size(FO.EventTimes,2) == 2 %Events with start/stops
    %Events that end after the window start or start before the window end
    inwinevents = FO.EventTimes(:,2) >=thiseventwin(1) & FO.EventTimes(:,1)<=thiseventwin(2);
    
    inwineventtimes = mean(FO.EventTimes(inwinevents,:),2); %Take the midpoint
    inwinstarts = FO.EventTimes(inwinevents,1);
    inwinstops = FO.EventTimes(inwinevents,2);

    plot(FO.viewwin,inwinstarts,zeros(size(inwinstarts))-500,...
        'o','color',[0 0.6 0],'markersize',5,'linewidth',3)
    plot(FO.viewwin,inwinstops,zeros(size(inwinstops))-500,...
        'o','color',[0.6 0 0],'markersize',5,'linewidth',3)
    
else
    error('EventTimes error, call in the big guns!')
end

%Plot the behavior
if isfield(FO,'behavior')
    numbeh = length(FO.behavior);
    for bb = 1:numbeh
        inwinbehavior = InIntervals(FO.behavior{bb}.timestamps,thiseventwin);
        behnorm = bz_NormToRange(FO.behavior{bb}.data,ywinrange);
        plot(FO.behavior{bb}.timestamps(inwinbehavior),behnorm(inwinbehavior))
    end
end

%Passthrough info from the plot
viewinfo.thiseventwin = thiseventwin;
viewinfo.inwinevents = inwineventtimes;

%UPdate event number display
set(FO.thiseventdisplay,'String',round(FO.currevent))

%% Deal with Flags and Comments
%Update Comment/Flag Window to reflect current event
%these try statements are to deal with FlagsAndComments not being made yet (do better)
if isfield(FO,'FlagsAndComments')
    if isfield(FO.FlagsAndComments,(FO.viewmode))
        %Has it already been flagged?
        [isflagged,flagidx] = ismember(FO.currevent,FO.FlagsAndComments.(FO.viewmode).flags);
        switch isflagged
            case true
                set(FO.eventcomment,'String',FO.FlagsAndComments.(FO.viewmode).comments{flagidx})
                set(FO.flageventbutton,'String','Unflag')
            case false
                set(FO.eventcomment,'String','Event Comments')
                set(FO.flageventbutton,'String','Flag')        
        end
        
    else %This is ugly :(
        set(FO.eventcomment,'String','Event Comments')
        set(FO.flageventbutton,'String','Flag')
    end
else %If no flags/comments have been created yet
    set(FO.eventcomment,'String','Event Comments')
    set(FO.flageventbutton,'String','Flag')
end
%% Focus on window
%This is a terrible fix to the figure focus problem that slows down window
%switching. sad. only necessary if focus has been moved away from current

% set(0,'currentfigure',FO.fig); %These is supposed to fix the post-button resize bug... 
% set(FO.fig,'currentaxes',FO.viewwin); %but do not. sad.
% figure(FO.fig)
% set(gcf,'CurrentObject',gcf)

%figure
set(findobj(FO.fig,'Type','uicontrol'),'Enable','off');
drawnow;
set(findobj(FO.fig,'Type','uicontrol'),'Enable','on');
end

