function metaInfoFigure(filebase)
% no inputs, assumes you're in directory with .dat, .meta, .tsp and .mpg files

%% get filebase if doesn't exist, save
if ~exist('filebase','var')
    choice = questdlg('No basename entered, use current directory name as file basename?','No Basename','Yes','Cancel','Yes');
    if strmatch(choice,'Cancel')
        return
    elseif strmatch(choice,'Yes')
         [pathstr, name, ext]=fileparts(cd);
         filebase = name;
    end
end

gd.filebase = filebase;


%% Get or create meta-info about the recording sessions in this folder
if FileExists([filebase, '-metaInfo.mat']);
    load([filebase, '-metaInfo.mat']);
else
    varargout = getDatAndVideoSessionInfo(filebase);
    metaInfo = varargout.metaInfo;
end


%% Make Figure
gd.metaInfoFig = figure('Position',[50 50 700 920],'Units','Pixels');

%% Setting up display of file durations
numSessions = size (metaInfo.names,2);

gd.displayRegionPanel = uipanel('Units','Pixels','Position',[7 95 686 838]);

topui = 890;

gd.columnTitlesNames = uicontrol('Style','text','String','File Name','Units','Pixels','Position',[15 topui 140 20]);
gd.columnTitlesStartSec = uicontrol('Style','text','String','Start(s)','Units','Pixels','Position',[160 topui 60 20]);
gd.columnTitlesStopSec = uicontrol('Style','text','String','Stop(s)','Units','Pixels','Position',[225 topui 60 20]);
gd.columnTitlesDurationSec = uicontrol('Style','text','String','Duration(s)','Units','Pixels','Position',[290 topui 70 20]);
gd.columnTitlesStartMinSec = uicontrol('Style','text','String','Start(m:s)','Units','Pixels','Position',[365 topui 70 20]);
gd.columnTitlesStopMinSec = uicontrol('Style','text','String','Stop(m:s)','Units','Pixels','Position',[440 topui 70 20]);
gd.columnTitlesDurationMinSec = uicontrol('Style','text','String','Duration(m:s)','Units','Pixels','Position',[515 topui 90 20]);

for a = 1:numSessions;
    gd.nameBoxes(a) = uicontrol('Style','text','String',metaInfo.names{a},'Units','Pixels','Position',[15 topui-(a*25) 140 20]);
    gd.startSec(a) = uicontrol('Style','text','String',num2str(metaInfo.startTimeInSec{a}),'Units','Pixels','Position',[160 topui-(a*25) 60 20]);
    gd.stopSec(a) = uicontrol('Style','text','String',num2str(metaInfo.stopTimeInSec{a}),'Units','Pixels','Position',[225 topui-(a*25) 60 20]);
    gd.durationSec(a) = uicontrol('Style','text','String',num2str(metaInfo.durationInSec{a}),'Units','Pixels','Position',[290 topui-(a*25) 70 20]);
    
    startmin = floor(metaInfo.startTimeInSec{a}/60);
    startsec = rem(metaInfo.startTimeInSec{a},60);
    startstring = [num2str(startmin),':',num2str(startsec)];
        
    stopmin = floor(metaInfo.stopTimeInSec{a}/60);
    stopsec = rem(metaInfo.stopTimeInSec{a},60);
    stopstring = [num2str(stopmin),':',num2str(stopsec)];
    
    durationmin = floor(metaInfo.durationInSec{a}/60);
    durationsec = rem(metaInfo.durationInSec{a},60);    
    durationstring = [num2str(durationmin),':',num2str(durationsec)];
    
    gd.startMinSec(a) = uicontrol('Style','text','String',startstring,'Units','Pixels','Position',[365 topui-(a*25) 70 20]);
    gd.stopMinSec(a) = uicontrol('Style','text','String',stopstring,'Units','Pixels','Position',[440 topui-(a*25) 70 20]);
    gd.durationMinSec(a) = uicontrol('Style','text','String',durationstring,'Units','Pixels','Position',[515 topui-(a*25) 90 20]);
end
    
%% Creating box to allow user to reference to movie timeframes
gd.inputRegionPanel = uipanel('Units','Pixels','Position',[7 3 686 89]);

gd.inputRegionTitle = uicontrol('Style','text','FontWeight','bold','String','FIND MOVIE TIMESTAMP FOR A GIVEN RECORDING FILE TIMESTAMP:','Units','Pixels','Position',[15 60 500 20]); 
gd.inputMinSecBoxTitle = uicontrol('Style','text','String','Enter time in min:sec','FontWeight','bold','Units','Pixels','Position',[10 35 160 20]); 
gd.inputSecBoxTitle = uicontrol('Style','text','String','Enter time in seconds','FontWeight','bold','Units','Pixels','Position',[10 10 160 20]); 

gd.inputMinSecEditBox = uicontrol('Style','edit','Units','Pixels','Position',[180 35 60 20],'Callback',@inputMinSecCallback); 
gd.inputSecEditBox = uicontrol('Style','edit','Units','Pixels','Position',[180 10 60 20],'Callback',@inputSecCallback); 

gd.outputMovieTitleLabel = uicontrol('Style','text','String','Moviename','Units','Pixels','Position',[250 35 130 20]); 
gd.outputMovieSecLabel = uicontrol('Style','text','String','Movie sec','Units','Pixels','Position',[400 35 70 20]); 
gd.outputMovieMinSecLabel = uicontrol('Style','text','String','Movie m:s','Units','Pixels','Position',[480 35 90 20]); 
gd.outputMovieTitleBox = uicontrol('Style','text','Units','Pixels','Position',[250 15 130 20]); 
gd.outputMovieSecBox = uicontrol('Style','text','Units','Pixels','Position',[400 15 70 20]); 
gd.outputMovieMinSecBox = uicontrol('Style','text','Units','Pixels','Position',[480 15 90 20]); 


gd.metaInfo = metaInfo;
guidata(gd.metaInfoFig,gd);

function inputMinSecCallback(obj,ev)
gd = guidata(obj);
minsec = get(gd.inputMinSecEditBox,'String');
errorstring = 'You must enter time in format min:sec (ie with colon in between)';
colonplace = findstr(':',minsec);
numchar = length(minsec);

if isempty(colonplace); %if no colon, assume all minutes
    minutes = str2num(minsec);
    seconds = 0;
    totalseconds = minutes*60+seconds;
else
    minutes = str2num(minsec(1:(colonplace-1)));
    seconds = str2num(minsec((colonplace+1):end));
    totalseconds = minutes*60+seconds;
end

[moviename, moviesec] = findMovieTimeStamp(gd.filebase,totalseconds,gd.metaInfo);
set(gd.outputMovieTitleBox,'String',moviename)
set(gd.outputMovieSecBox,'String',num2str(moviesec))
moviemin = floor(moviesec/60);
moviesec = rem(moviesec,60);
movieminsecstring = [num2str(moviemin),':',num2str(moviesec)];
set(gd.outputMovieMinSecBox,'String',movieminsecstring);
set(gd.inputSecEditBox,'String',[]); 

%%
function inputSecCallback(obj,ev)
gd = guidata(obj);
seconds = get(gd.inputSecEditBox,'String');
totalseconds = str2num(seconds);

[moviename, moviesec] = findMovieTimeStamp(gd.filebase,totalseconds,gd.metaInfo);
set(gd.outputMovieTitleBox,'String',moviename)
set(gd.outputMovieSecBox,'String',num2str(moviesec))
moviemin = floor(moviesec/60);
moviesec = rem(moviesec,60);
movieminsecstring = [num2str(moviemin),':',num2str(moviesec)];
set(gd.outputMovieMinSecBox,'String',movieminsecstring);
set(gd.inputMinSecEditBox,'String',[]); 

%%
function [moviename, moviesec] = findMovieTimeStamp(filebase,sec,metaInfo)

for a = 1:size(metaInfo.names,2)
    afterstart(a) = metaInfo.startTimeInSec{a}<sec;
    beforestop(a) = metaInfo.stopTimeInSec{a}>sec;
end
movieindex = find(afterstart.*beforestop);
moviename = metaInfo.names{movieindex};

if FileExists([moviename, '-tsp.mat']);
    tsp = load([moviename, '-tsp.mat']);
    tsp = tsp.tsp;
else
    tsp = load([moviename, '.tsp']);
    save([cd,'/', moviename, '-tsp.mat'], 'tsp');
    ['Saved TSP mat file'];
end
tsp = tsp(:,1); %for now, no need to keep LED info.

tstampDur = metaInfo.stopTimeStampInMs{movieindex} - metaInfo.startTimeStampInMs{movieindex};
datDur = double(metaInfo.fileBytes{movieindex}/(metaInfo.numChannels{movieindex}*2*20));

tsp = tsp - metaInfo.startTimeStampInMs{movieindex};
tsp = tsp*(datDur/tstampDur);
tsp = tsp/1000;
% f = sec:(1/30):(sec +  duration);

thisfilesec = sec - metaInfo.startTimeInSec{movieindex};

framenum = dsearchn(tsp, thisfilesec);
moviesec = tsp(framenum);