function [ DetectionReview ] = DetectionReview(obj,event )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%obj = findobj('tag','EventExplorerMaster');  
FO = guidata(obj); 
%% Select the time windows to look at 
%User input for this
numwins = 30; %number of windows to look at. determine to maximize sampling or have user input (FO.uinput.field?)
%Selecting from events:
%randevents = randsample(FO.EventTimes,numevents);
%Selecting from random times (RestrictInts takes way too long...)
set(findobj(FO.fig,'Type','uicontrol'),'Enable','off');
drawnow;
if ~isempty(FO.detectionints)
    [status,interval,index] = InIntervals(FO.data.lfp.timestamps,double(FO.detectionints));
    restrictedtimes = FO.data.lfp.timestamps(status);
else
    restrictedtimes = FO.data.lfp.timestamps;
end

%use InIntervals
set(findobj(FO.fig,'Type','uicontrol'),'Enable','on');
randevents = randsample(restrictedtimes,numwins);
%Find any events that are within winsize of another event to remove them?
%Also should look at windows in which no events were detected?  Maybe just
%random windows in the detectionwin
closeevents = abs(diff(randevents))<FO.winsize;
% numtries =0;
% while any(closeevents) && numtries<50
%     randevents(closeevents)
% end
%Problem:overlapping ints
%[ ints ] = MergeSeparatedInts( ints,minseparation )


%%
%Store the current user action
FO.currentuseraction = 'MarkEvents';
FO.markedevents = []; %clear any previously marked events
FO.viewmode = 'timepoint';
guidata(FO.fig, FO); %store the data in the figure

counter = numwins;
miss=[];hit=[];falsealarm=[];
lookedatwins = [];

%UI panel for event review
FO.EventPanel = uipanel('Title','Detection Review','FontSize',12,...
        'Position',[.65 .05 0.25 0.32]);
%Buttons for next loop or finishing
nextbtn = uicontrol('Parent',FO.EventPanel,...
    'Position',[200 20 150 40],'String','Next Window (Return)',...
          'Callback','uiresume(gcbf)');
qutbtn = uicontrol('Parent',FO.EventPanel,...
    'Position',[200 65 150 40],'String','Quit Early (Esc)',...
         'Callback','QUITLOOP=true;uiresume(gcbf)');
%Instruction text
instructtext = uicontrol('Parent',FO.EventPanel,...
    'Position',[25 100 400 100],'style','text',...
    'string',{'Instructions:',...
    '   - LEFT click missed events (o)','   - RIGHT click False Alarms (x)'},...
    'HorizontalAlignment','left'); 
%Display counter text
correcttext = uicontrol('Parent',FO.EventPanel,...
    'Position',[25 60 125 15],'style','text',...
    'string','Correct: 0','HorizontalAlignment','left'); 
misstext = uicontrol('Parent',FO.EventPanel,...
    'Position',[25 40 125 15],'style','text',...
    'string','Miss: 0','HorizontalAlignment','left');
FAtext = uicontrol('Parent',FO.EventPanel,...
    'Position',[25 20 125 15],'style','text',...
    'string','FA: 0','HorizontalAlignment','left');
countetext = uicontrol('Parent',FO.EventPanel,...
    'Position',[25 100 175 15],'style','text',...
    'string',['Windows Remaining: ',num2str(counter)],'HorizontalAlignment','left');

%The Event Review Loop
QUITLOOP = false;
while counter>0 && QUITLOOP~=true
    %EventUI(FO)   %function that will run the user interface
    thiseventtime = randevents(counter);  
    FO.currevent = thiseventtime; guidata(FO.fig, FO);
    viewinfo = EventVewPlot;
    uiwait(FO.fig)  %Wait here until the user clicks quit or next
    
    lookedatwins = [lookedatwins; viewinfo.thiseventwin];
    inwinevents = viewinfo.inwinevents;
    
    
    %Tally the user selections for that window
    obj = findobj('tag','EventExplorerMaster');  FO = guidata(obj);
    if ~isempty(FO.markedevents)
        FO.markedevents(isnan(FO.markedevents(:,1)))=[]; %remove nans - bug fix
        miss = [miss; FO.markedevents(FO.markedevents(:,3)==1,1)];
        %This is messy to account for interp errors with 0-1 reference points
        if isempty(FO.markedevents(FO.markedevents(:,3)==3,1))
            newfalsealarms = [];
        elseif length(FO.markedevents(FO.markedevents(:,3)==3,1))==1 && length(inwinevents)==1
            newfalsealarms = inwinevents;
        else
            newfalsealarms = interp1(inwinevents,inwinevents,...
                FO.markedevents(FO.markedevents(:,3)==3,1),'nearest');
        end
        falsealarm = [falsealarm; newfalsealarms];
    end
    hit = [hit; inwinevents(~ismember(inwinevents,falsealarm))];
    FO.markedevents = [];  guidata(FO.fig, FO); %reset the marked events
    
    counter = counter-1;
    set(correcttext, 'String',['Correct: ',num2str(length(hit))]);
    set(misstext, 'String', ['Miss: ',num2str(length(miss))]);
    set(FAtext, 'String', ['FA: ',num2str(length(falsealarm))]);
    set(countetext, 'String', ['Windows Remaining: ',num2str(counter)]);
end
%Calculate total number of miss,hit,FA
numMiss = length(miss);
numHit = length(hit);
numFA = length(falsealarm);    

%Calculate total amount of time/percentage of detection time (detectionintervals) looked at

%Put things in the output structure
DetectionReview.lookedatwins = lookedatwins;
DetectionReview.miss = miss; 
DetectionReview.hit = hit;
DetectionReview.falsealarm = falsealarm;
DetectionReview.estMissperc = numMiss./(numHit+numMiss);
DetectionReview.estFAperc = numFA./(numHit+numFA);
DetectionReview.ReviewDate = today;
DetectionReview.EventsType = FO.EventName;

%UI: Done!  Would you like to save the results to (eventsfilename?)
%Make function that does this: SaveResults(FO,EEoutput)
if isfield(FO,'eventsfilename')
    button = questdlg(['Event Review Complete! Would you like to add the results to ',...
        FO.eventsfilename,'?'],'Good Job!');
    switch button
        case 'Yes'
            %Load the events file, add the field, save the events file
            try %Only do this if the correct named structure lives in the file
                eventsfile = load(FO.eventsfilename,FO.EventName);
                eventsfile.(FO.EventName).EventExplorer.DetectionReview = DetectionReview;
                save(FO.eventsfilename,'-struct','eventsfile',FO.EventName,'-append')
            catch
                warndlg({' Save failed... ',[FO.eventsfilename,' may not ',...
                    'contain a structure titled ',FO.EventName,'.'],...
                    'Or you may not have sudo priviliges...?'},'Oh No!')
            end
    end
end

%% Return to explorer mode
FO.viewmode = 'events';
FO.currentuseraction = 'none';
FO.currevent = 1;
set(FO.EventPanel,'Visible','off')
FO.DetectionReview = DetectionReview;
guidata(FO.fig, FO); %Save the detection review back to GUI data
end

