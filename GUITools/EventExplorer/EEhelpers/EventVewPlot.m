function [ viewinfo ] = EventVewPlot
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
obj = findobj('tag','EventExplorerMaster');  FO = guidata(obj); 

%Two options of current timepoint: Event Number, arbitraty timepoint
switch FO.viewmode
    case 'timepoint' 
        timepoint = FO.currevent;
    case 'events'
        timepoint = FO.EventTimes(FO.currevent);
    case 'FAs'
        timepoint = FO.DetectionReview.falsealarm(FO.currevent);
    case 'misses'
        timepoint = FO.DetectionReview.miss(FO.currevent);
end

%Get the view window and events inside it
thiseventwin = timepoint+FO.winsize.*[-0.5 0.5];
inwinevents = FO.EventTimes(FO.EventTimes>=thiseventwin(1) &FO.EventTimes<=thiseventwin(2));

%Plot
%viewwin = subplot(3,1,2,'ButtonDownFcn',@MouseClick);
%set(gca,'ButtonDownFcn', @MouseClick)
hold(FO.viewwin,'off')
bz_MultiLFPPlot(FO.data.lfp,'timewin',thiseventwin,'spikes',FO.data.spikes,...
    'axhandle',FO.viewwin,'scaleLFP',FO.scaleLFP)
hold on
plot(FO.viewwin,inwinevents,zeros(size(inwinevents)),'o','color',[0 0.6 0])


%Passthrough info from the plot
viewinfo.thiseventwin = thiseventwin;
viewinfo.inwinevents = inwinevents;

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

