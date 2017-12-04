function [SleepState,durationparams] = StatesToFinalScoring(WAKEints,NREMints,REMints)
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

[episodeintervals{1}] = IDStateEpisode(WAKEints,maxWAKEEpisodeInterruption,minWAKEEpisodeDuration);
[episodeintervals{2}] = IDStateEpisode(NREMints,maxNREMEpisodeInterruption,minNREMEpisodeDuration);
[episodeintervals{3}] = IDStateEpisode(REMints,maxREMEpisodeInterruption,minREMEpisodeDuration);

%% Exclude REM from NREM episodes that had bridged non-NREM periods!
OrigNREMEpisodes = episodeintervals{2};%this will be a list of acceptable episodes
AllNewNREMEpisodes = [];%will populate this list with snippets of interrupted episodds

REMints = sortrows(REMints);
[~,NREMwRemInside] = InIntervals(REMints(:,1),episodeintervals{2});

NREMwRemInside(NREMwRemInside == 0) = [];
badNrem = unique(NREMwRemInside);
for nix = length(badNrem):-1:1%backwards bc will be cutting out episodes
    tnrem = badNrem(nix);
    thisNREMint = episodeintervals{2}(tnrem,:);
    OrigNREMEpisodes(tnrem,:) = [];%delete bad episodes from permanent list
    
    whichREMinThisNREM = InIntervals(REMints(:,1),thisNREMint);
    whichREMinThisNREM = find(whichREMinThisNREM);
    for rix = 1:length(whichREMinThisNREM)
        tREM = whichREMinThisNREM(rix);
        newNREMint = [thisNREMint(1,1) REMints(tREM,1)];%keep segment of NREM before start of REM
        AllNewNREMEpisodes = cat(1,AllNewNREMEpisodes,newNREMint);
        
        thisNREMint = [REMints(tREM,2) thisNREMint(1,2)];%NREM remaining after REM, for next cycle
    end
    AllNewNREMEpisodes = cat(1,AllNewNREMEpisodes,thisNREMint);%toss on whatever's left at the end
end
NREMEpisodes = cat(1,OrigNREMEpisodes,AllNewNREMEpisodes);%combine new and old
NREMEpisodes = sortrows(NREMEpisodes);%sort by start time
NREMEpisodes(diff(NREMEpisodes,[],2)<minNREMEpisodeDuration,:) = [];%get rid of too-short episodes
episodeintervals{2} = NREMEpisodes;

%% Identify MAs preceeding REM and separate them into MA_REM
%Find the MA states that are before (or between two) REM states
[ ~,~,MAInREM,~ ] = FindIntsNextToInts(MAIntervals,REMints);
%"Real" MAs are those that are not before REM state
realMA = setdiff(1:size(MAIntervals,1),MAInREM);

MAInREM = MAIntervals(MAInREM,:);
MAIntervals = MAIntervals(realMA,:);

%% Save: buzcode format
SleepState.intsRaw.NREMstate = NREMints;
SleepState.intsRaw.REMstate = REMints;
SleepState.intsRaw.WAKEstate = WAKEIntervals;

SleepState.intsWatson2016.WAKEpisodes = episodeintervals{1};
SleepState.intsWatson2016.NREMepisodes = episodeintervals{2};
SleepState.intsWatson2016.REMepisodes = episodeintervals{3};
SleepState.intsWatson2016.NREMpackets = packetintervals;
SleepState.intsWatson2016.MAstate = MAIntervals;
SleepState.intsWatson2016.MAInREM = MAInREM;





