function varargout = getDatAndVideoSessionInfo(fb,forceRegenerate)

% TSPS = dir('*tsp');
metas = dir('*.meta');
if ~exist('forceRegenerate','var')
    forceRegenerate = 0;
end

if length(metas) > 1
    if forceRegenerate
        if FileExists([fb, '-metaInfo.mat'])
            movefile([fb, '-metaInfo.mat'],[fb, '-metaInfo_Pre' date '.mat'])
        end
    end
    if ~FileExists([fb, '-metaInfo.mat'])
        makeOrderDiag(metas, fb);
    end
    metaInfo = load([fb, '-metaInfo.mat']);
%     t = s.sessOrder.times;
%     
%     rSess = [];
%     for i = 1:size(t, 1)
%         if time >= t(i, 1) & time <= t(i, 2)
%             rSess = i;
%         end
%     end
%     if isempty(rSess)
%         error('Requested time not found, the managment regets any inconvenience');
%     end
%     name1  = s.sessOrder.names{rSess};
%     time1 = time - t(rSess, 1);
%     offset1 = t(rSess, 1);
% else
%     name1 = TSPS(1).name(1:end - 4);
%     time1 = time;
%     offset1 = 0;
    if nargout>0
        varargout{1} = metaInfo;
    end
else
    
end
% 
% playVideoSyncChunkIN(time1, duration, name1, offset1);
end

function makeOrderDiag(metas, fb)
fig1 = figure;
lists1 = {};
choices = [];
for i = 1:length(metas)
    choices = [choices, int2str(i), '|'];
end
choices = choices(1:(end - 1));

for i = 1:length(metas)
    lists1{i} = uicontrol('Style', 'popup', 'String', choices, 'Position', [ 20, 400 - (i*40),  60,  20],'Value',i); %defaults to alphabetical (sometimes same as numerical) order
    uicontrol('Style', 'text', 'String', metas(i).name, 'HorizontalAlignment', 'left', 'Position', [ 90, 400 - (i*40), 300, 20]);
end

cbut = uicontrol('Style', 'pushbutton', 'String', 'CONFIRM', 'Position', [ 300    20    60    20], 'Callback', {@confirmMe, metas, lists1, fig1, fb});

uiwait(fig1);
close(fig1);
end


function confirmMe(x, y, metas, lists1, fig1, fb)

v1 = [];
c = [];
for i = 1:length(lists1)
    v1 = [v1, get(lists1{i}, 'Value')];
    c = [c, i];
end

if sum(sort(v1) == c) == length(v1)
    metas = metas(v1);
    metaInfo = [];
%     sessMetaData.names = {};
%     metaInfo.startStopTimes = [];
    last = 0;
    for i = 1:length(v1)
        n1 = [];
        ChanNum = [];
        SamplesPerSec = [];
        DatSize = [];
        StartTimestamp = [];
        EndTimestamp = [];
        startTimeByClock = [];
        startDate = [];

%         n1 = metas(i).name(1:(end - 4));%for .tsp
        n1 = metas(i).name(1:(end - 5));%for .meta
        if i == 9;
            1;
        end
        fid=fopen([n1, '.meta']);
        tline= fgetl(fid);
        while ischar(tline)
            try
                if strcmp(tline(1:20),'TimeStamp of the end')
                    tline=tline(59:end);
                    EndTimestamp=sscanf(tline,'%s',1);
                    EndTimestamp=str2num(EndTimestamp);
                end
            catch end
            try
                if strcmp(tline(1:22),'TimeStamp of the start')
                    tline=tline(61:end);
                    StartTimestamp=sscanf(tline,'%s',1);
                    StartTimestamp=str2num(StartTimestamp);
                end
            catch
            end
            try
                if strcmp(tline(1:9),'Number of')
                    tline=tline(31:end);
                    ChanNum=sscanf(tline,'%s',1);
                    ChanNum=str2num(ChanNum);
                end
            catch end
            try
                if strcmp(tline(1:17),'File size (bytes)')
                    tline=tline(21:end);
                    DatSize=sscanf(tline,'%s',1);
                    DatSize=str2num(DatSize);
                end
            catch end
            try
                if strcmp(tline(1:15),'Sampling rate =')
                    tline=tline(17:end);
                    SamplesPerSec=sscanf(tline,'%s',1);
                    SamplesPerSec=str2num(SamplesPerSec);
                end
            catch end
            try 
                if strcmp(tline(1:20),'Recording start time')
                    tline = tline(24:31);
                    startTimeByClock = sscanf(tline,'%s',1);
                end
            catch end
            try 
                if strcmp(tline(1:20),'Recording start date')
                    tline = tline(28:38);
                    startDate = tline;
                end
            catch end

            tline= fgetl(fid);
        end
        fclose(fid);
        metaInfo.names{i} = n1;
        metaInfo.numChannels{i} = ChanNum;
        sessMetadata.samplingRate{i} = SamplesPerSec;
        metaInfo.fileBytes{i} = DatSize;
        metaInfo.startTimeStampInMs{i} = StartTimestamp;
        metaInfo.stopTimeStampInMs{i} = EndTimestamp;
        metaInfo.startTimeByClock{i} = startTimeByClock;
        metaInfo.startDate{i} = startDate;
        
        d = double(DatSize/(2*ChanNum*SamplesPerSec));
        metaInfo.startTimeInSec{i} = last;
        last = last + d;
        metaInfo.stopTimeInSec{i} = last;
        metaInfo.durationInSec{i} = d;
    end
    %    sessOrder.times = (sessOrder.times - min(min(sessOrder.times)))/1000;
    save([cd,'/', fb, '-metaInfo.mat'], 'metaInfo');
    uiresume(fig1);
else
    disp('Error: Ordering does not appear to make sense')
end
end



