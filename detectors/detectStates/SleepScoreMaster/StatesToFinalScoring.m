function [SleepState,durationparams] = StatesToFinalScoring(NREMints,WAKEints,REMints)

minPacketDuration = 30;
maxMicroarousalDuration = 100;

maxEpisodeDuration = 40;
minNREMEpisodeDuration = 20;
minWAKEEpisodeDuration = 20;
minREMEpisodeDuration = 20;

durationparams.minPacketDuration = minPacketDuration;
durationparams.maxMicroarousalDuration = maxMicroarousalDuration;
durationparams.maxEpisodeDuration = maxEpisodeDuration;
durationparams.minNREMEpisodeDuration = minNREMEpisodeDuration;
durationparams.minWAKEEpisodeDuration = minWAKEEpisodeDuration;
durationparams.minREMEpisodeDuration = minREMEpisodeDuration;


NREMlengths = NREMints(:,2)-NREMints(:,1);
packetintervals = NREMints(NREMlengths>=minPacketDuration,:);

WAKElengths = WAKEints(:,2)-WAKEints(:,1);
MAIntervals = WAKEints(WAKElengths<=maxMicroarousalDuration,:);
WAKEIntervals = WAKEints(WAKElengths>maxMicroarousalDuration,:);

[episodeintervals{2}] = IDStateEpisode(NREMints,maxEpisodeDuration,minNREMEpisodeDuration);
[episodeintervals{1}] = IDStateEpisode(WAKEints,maxEpisodeDuration,minWAKEEpisodeDuration);
[episodeintervals{3}] = IDStateEpisode(REMints,maxEpisodeDuration,minREMEpisodeDuration);


%% Identify MAs preceeding REM and separate them into MA_REM
%Find the MA states that are before (or between two) REM states
[ ~,~,MA_REM,~ ] = FindIntsNextToInts(MAIntervals,REMints);
%"Real" MAs are those that are not before REM state
realMA = setdiff(1:size(MAIntervals,1),MA_REM);

MA_REM = MAIntervals(MA_REM,:);
MAIntervals = MAIntervals(realMA,:);



% %% Save: old format - remove once compatible with StateEditor
% StateIntervals.NREMstate = NREMints;
% StateIntervals.REMstate = REMints;
% StateIntervals.WAKEstate = WAKEIntervals;
% StateIntervals.NREMepisode = episodeintervals{2};
% StateIntervals.REMepisode = episodeintervals{3};
% StateIntervals.WAKEepisode = episodeintervals{1};
% StateIntervals.NREMpacket = packetintervals;
% StateIntervals.MAstate = MAIntervals;
% StateIntervals.MA_REM = MA_REM;

%% Save: buzcode format
SleepState.ints.NREMstate = NREMints;
SleepState.ints.REMstate = REMints;
SleepState.ints.WAKEstate = WAKEIntervals;
SleepState.ints.NREMepisode = episodeintervals{2};
SleepState.ints.REMepisode = episodeintervals{3};
SleepState.ints.WAKEepisode = episodeintervals{1};
SleepState.ints.NREMpacket = packetintervals;
SleepState.ints.MAstate = MAIntervals;
SleepState.ints.MA_REM = MA_REM;





