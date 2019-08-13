function [zeitgeberTime relativeTime absTime] = bz_getZeitgeberTime(basename,downSampleFactor)
% assumes pwd is basepath for now


%% get timestamps for recording starts
parts = strsplit(basename, '_');
date = parts{end-1};

if strcmp(parts{end},'merge') % if multiple sessions are merged, find the start times of each
    fileList = dir('*');
    
    for f = 1:length(fileList)
        if fileList(f).isdir
        for p = 1:length(parts)
            if ~isempty(strfind(fileList(f).name,parts{p}))
               matches(f,p) = 1;
            else
               matches(f,p) = 0;
            end
        end
        end
    end    
    idx = find(sum(matches')==length(parts)-1);  % exlude '_merge' string
    for ind = 1:length(idx)
        temp = strsplit(fileList(idx(ind)).name,'_');
        recordingStart{ind} = temp{end};
    end
else
    recordingStart{1} = parts{end};
end

for ind = 1:length(recordingStart)
timeSeconds(ind) = int32(str2num(recordingStart{ind}(1:2)) * 60 * 60 + ...
                   str2num(recordingStart{ind}(3:4)) * 60 + ...
                   str2num(recordingStart{ind}(5:6)));
end

%% now get relative time from time.dat file
fileinfo = dir('time.dat');
num_samples = fileinfo.bytes/4; % int32 = 4 bytes
fid = fopen('time.dat', 'r');
t = fread(fid, num_samples, 'int32=>int32')';
t = downsample(t,downSampleFactor);
fclose(fid);

[a b] = min(t);  %% fixes output for time.dat files that 'wrap' at the int32 limit
if a < 0
    t = double(t)./20000;
    for i=b:length(t)
    t(i) = t(i-1)+1;
    end
else
    t = single(t)./20000;
end
% t = t / 20000;

idx = find(t==0);
idx(length(idx)+1) = length(t); % add end ts

for ind = 1:length(idx)-1
    zeitgeberTime(idx(ind):idx(ind+1)) = single(t(idx(ind):idx(ind+1)))/20000 + single(timeSeconds(ind));
end

zeitgeberTime = single(wrap(double(zeitgeberTime)/86400*pi*2,2))/ (2*pi) * 86400; % wrap it up

relativeTime = t;  % recording time

% zeitgeberTime = downsample(zeitgeberTime,downSampleFactor);
% relativeTime = downsample(relativeTime,downSampleFactor

% below is error prone if gaps in relative time, but WAY faster
absTime = linspace(datenum(2000 + str2num(date(1:2)),str2num(date(3:4)),str2num(date(5:6)),...
        str2num(parts{end}(1:2)),str2num(parts{end}(3:4)),double(relativeTime(1))),...
        datenum(2000 + str2num(date(1:2)),str2num(date(3:4)),str2num(date(5:6)),...
        str2num(parts{end}(1:2)),str2num(parts{end}(3:4)),double(relativeTime(end))),length(relativeTime));
    
% switch to this for accuracy
% for ind = 1:length(relativeTime)
%     absTime(ind) = datenum(2000 + str2num(date(1:2)),str2num(date(3:4)),str2num(date(5:6)),...
%         str2num(parts{end}(1:2)),str2num(parts{end}(3:4)),double(relativeTime(ind)));
% end
























