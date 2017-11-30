function [SleepState,durationparams] = StatesToFinalScoring(NREMints,WAKEints,REMints)
% Takes sleep state data as a series of raw intervals and imposes max and
% min durations and interruption criteria to make WAKE, NREM, REM Episodes,
% Packets, Microarousals and Microarusals in REM.
%
% INPUTS
%
% OUTPUTS
%
% Called in TheStateEditor saveStates function
% Dan Levenstein & Brendon Watson 2016

minPacketDuration = 30;
minWAKEEpisodeDuration = 20;
minNREMEpisodeDuration = 20;
minREMEpisodeDuration = 20;

maxMicroarousalDuration = 100;

maxWAKEEpisodeInterruption = 40;
maxNREMEpisodeInterruption = maxMicroarousalDuration;
maxREMEpisodeInterruption = 40;

durationparams = v2struct(minPacketDuration, minWAKEEpisodeDuration,...
    minNREMEpisodeDuration, minREMEpisodeDuration,...
    maxMicroarousalDuration,...
    maxWAKEEpisodeInterruption, maxNREMEpisodeInterruption,...
    maxREMEpisodeInterruption);

NREMlengths = NREMints(:,2)-NREMints(:,1);
packetintervals = NREMints(NREMlengths>=minPacketDuration,:);

WAKElengths = WAKEints(:,2)-WAKEints(:,1);
MAIntervals = WAKEints(WAKElengths<=maxMicroarousalDuration,:);
WAKEIntervals = WAKEints(WAKElengths>maxMicroarousalDuration,:);

[episodeintervals{2}] = IDStateEpisode(NREMints,maxNREMEpisodeInterruption,minNREMEpisodeDuration);
[episodeintervals{1}] = IDStateEpisode(WAKEints,maxWAKEEpisodeInterruption,minWAKEEpisodeDuration);
[episodeintervals{3}] = IDStateEpisode(REMints,maxREMEpisodeInterruption,minREMEpisodeDuration);

%% Identify MAs preceeding REM and separate them into MA_REM
%Find the MA states that are before (or between two) REM states
[ ~,~,MA_REM,~ ] = FindIntsNextToInts(MAIntervals,REMints);
%"Real" MAs are those that are not before REM state
realMA = setdiff(1:size(MAIntervals,1),MA_REM);

MA_REM = MAIntervals(MA_REM,:);
MAIntervals = MAIntervals(realMA,:);

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





