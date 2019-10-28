function RecordingSecondsToTimeSeconds(basepath)
% Store correspondences between recoring (ie dat file) second timestamps
% and clock time.  0.1sec resolution. 
% Can only be run after bz_DatFileMetadata.m and TimeFromLightCycleStart.m
% Brendon Watson 2016

if ~exist('basepath','var')
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);

%% Check if already bins made
if exist(fullfile(basepath,[basename '_RecordingSecondVectors.mat']),'file')
    disp([fullfile(basepath,[basename '_RecordingSecondVectors.mat']) ' already exists, not re-making'])
end


%% Make bins if not already there
% load(fullfile(basepath,[basename '_DatInfo.mat']))
load(fullfile(basepath,[basename '_DatsMetadata.mat']))
load(fullfile(basepath,[basename '_SecondsFromLightsOn.mat']))

RecordingSeconds = .1:.1:sum(DatsMetadata.Recordings.Seconds);

% Clock seconds since light on (on day 1)
FromLightOn_ByClockSeconds = [];
RecordingStartsFromLightOnByClock = SecondsAfterLightCycleStart_PerFile;
RecordingEndsFromLightOnByClock = SecondsAfterLightCycleStart_PerFile+DatsMetadata.Recordings.Seconds;
for a = 1:length(RecordingStartsFromLightOnByClock)
    FromLightOn_ByClockSeconds = cat(2,FromLightOn_ByClockSeconds,RecordingStartsFromLightOnByClock(a)+.1:.1:RecordingEndsFromLightOnByClock(a));
end

% Clock seconds since start of file
FromRecordingOn_ByClockSeconds = FromLightOn_ByClockSeconds-SecondsAfterLightCycleStart_PerFile(1);
RecordingStartsFromRecordingOnByClock = RecordingStartsFromLightOnByClock-SecondsAfterLightCycleStart_PerFile(1);
RecordingEndsFromRecordingOnByClock = RecordingEndsFromLightOnByClock-SecondsAfterLightCycleStart_PerFile(1);


RecordingSecondVectors = v2struct(RecordingSeconds,...
    FromRecordingOn_ByClockSeconds,...
    FromLightOn_ByClockSeconds,...
    RecordingStartsFromLightOnByClock,RecordingEndsFromLightOnByClock,...
    RecordingStartsFromRecordingOnByClock,RecordingEndsFromRecordingOnByClock);

save(fullfile(basepath,[basename '_RecordingSecondVectors.mat']),'RecordingSecondVectors')

