function [SleepState,durationparams] = StatesToFinalScoring(NREMints,WAKEints,REMints, INTERints)

minPacketDuration = 30;
maxMicroarousalDuration = 100;

maxEpisodeDuration = 40;
minSWSEpisodeDuration = 20;
minWAKEEpisodeDuration = 20;
minREMEpisodeDuration = 20;
minINTEREpisodeDuration = 20; %stephanie added this

durationparams.minPacketDuration = minPacketDuration;
durationparams.maxMicroarousalDuration = maxMicroarousalDuration;
durationparams.maxEpisodeDuration = maxEpisodeDuration;
durationparams.minSWSEpisodeDuration = minSWSEpisodeDuration;
durationparams.minWAKEEpisodeDuration = minWAKEEpisodeDuration;
durationparams.minREMEpisodeDuration = minREMEpisodeDuration;
durationparams.minINTEREpisodeDuration = minINTEREpisodeDuration; %added by stephanie


SWSlengths = NREMints(:,2)-NREMints(:,1);
packetintervals = NREMints(SWSlengths>=minPacketDuration,:);

WAKElengths = WAKEints(:,2)-WAKEints(:,1);
MAIntervals = WAKEints(WAKElengths<=maxMicroarousalDuration,:);
WAKEIntervals = WAKEints(WAKElengths>maxMicroarousalDuration,:);
INTERlengths = INTERints(:,2)-INTERints(:,1); %added by stephanie
INTERIntervals = INTERints(INTERlengths>maxMicroarousalDuration,:); %added by stephanie

[episodeintervals{2}] = IDStateEpisode(NREMints,maxEpisodeDuration,minSWSEpisodeDuration);
[episodeintervals{1}] = IDStateEpisode(WAKEints,maxEpisodeDuration,minWAKEEpisodeDuration);
[episodeintervals{3}] = IDStateEpisode(REMints,maxEpisodeDuration,minREMEpisodeDuration);
[episodeintervals{4}] = IDStateEpisode(INTERints,maxEpisodeDuration,minINTEREpisodeDuration); %stephanie added this


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
SleepState.ints.INTERstate = INTERIntervals;
SleepState.ints.NREMepisode = episodeintervals{2};
SleepState.ints.REMepisode = episodeintervals{3};
SleepState.ints.WAKEepisode = episodeintervals{1};
SleepState.ints.INTERepisode = episodeintervals{4};
SleepState.ints.NREMpacket = packetintervals;
SleepState.ints.MAstate = MAIntervals;
SleepState.ints.MA_REM = MA_REM;





