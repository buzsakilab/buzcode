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
        timepoint = FO.EventReview.falsealarm(FO.currevent);
    case 'misses'
        timepoint = FO.EventReview.miss(FO.currevent);
end


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
viewinfo.inwinevents = inwinevents;
viewinfo.thiseventwin = thiseventwin;


%Update Comment/Flag Window to reflect current event
%these try statements are to deal with FlagsAndComments not being made yet (do better)
try iscommented = ~isempty(FO.FlagsAndComments.(FO.viewmode).comments{FO.currevent});
catch; iscommented=false; end
if iscommented 
    set(FO.eventcomment,'String',FO.FlagsAndComments.(FO.viewmode).comments{FO.currevent})
else set(FO.eventcomment,'String','Event Comments')
end

try isflagged = ismember(FO.currevent,FO.FlagsAndComments.(FO.viewmode).flags);
catch; isflagged = false; end
if isflagged
    set(FO.flageventbutton,'String','Unflag')
else set(FO.flageventbutton,'String','Flag')
end

% set(0,'currentfigure',FO.fig); %These is supposed to fix the post-button resize bug... 
% set(FO.fig,'currentaxes',FO.viewwin); %but do not. sad.
% figure(FO.fig)
% set(gcf,'CurrentObject',gcf)
%This is a terrible fix to the figure focus problem that slows down window
%switching. sad. only necessary if focus has been moved away from current
%figure
set(findobj(FO.fig,'Type','uicontrol'),'Enable','off');
drawnow;
set(findobj(FO.fig,'Type','uicontrol'),'Enable','on');
end

