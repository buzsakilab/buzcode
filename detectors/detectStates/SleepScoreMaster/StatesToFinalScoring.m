function [SleepState,SleepStateEpisodes] = StatesToFinalScoring(NREMints,WAKEints,REMints,DetectorName,DetectionParams)
% Takes sleep state data as a series of raw intervals and imposes max and
% min durations and interruption criteria to make WAKE, NREM, REM Episodes,
% Packets, Microarousals and Microarusals in REM.
%
% INPUTS
% - WAKEints - startstop pairs (nx2 array) of second numbers of starts and
%   ends respectively of periods of wake-like state.  These can be broken
%   later into wake, microarousal and other mini wake-like states
% - NREMints - startstop pairs (nx2 array) of second numbers of starts and
%   ends respectively of periods of NREM state
% - REMints - startstop pairs (nx2 array) of second numbers of starts and
%   ends respectively of periods of REM state
% - DetectorName - function that detected NREM, WAKE, REM states... 
%   essentially the calling function in most cases, such as 
%   SleepScoreMaster or TheStateEditor
% - DetectionParams - Parameters used by  
% 
% OUTPUTS
% - SleepState - Struct array of relatively less processed epochs relative
%   to SleepStateEpisodes.  
%   - .ints contains start-stop pair intervals in subfields: WAKEstate, 
%         NREMstate and REMstate.
%       - NREMstate: NREMstate same as the input NREMints  
%       - REMstate: REMstate same as the input REMints  
%       - WAKEstate: input WAKEints that are longer than the microarousal duration (100s at the writing of 
%           this comment)
%
% 
% - SleepStateEpisodes - Struct array of relatively more processed epochs
%   with minimum durations and maximum interruptions imposed.  Fields below
%   - WAKEpisodes: WAKEints inputs with at least [minWAKEEpisodeDuration] 
%   long and without interruptions more than [maxWAKEEpisodeInterruption]
%   - NREMepisodes
%   - REMepisodes
%   - NREMpackets
%   - MAstate
%   - MAInREM 
% 

%
% Called in [SleepScoreMaster] and [TheStateEditor "saveStates" & "ReClusterStates_In" functions]  
% Dan Levenstein & Brendon Watson 2016

if ~exist('DetectorName','var')
    DetectorName = 'NotSpecified';
end
if ~exist('DetectionParams','var')
    DetectionParams = [];
end
minPacketDuration = 30;
minWAKEEpisodeDuration = 20;
minNREMEpisodeDuration = 20;
minREMEpisodeDuration = 20;

maxMicroarousalDuration = 100;

maxWAKEEpisodeInterruption = 40;
maxNREMEpisodeInterruption = maxMicroarousalDuration;
maxREMEpisodeInterruption = 40;

detectionparams_episodes = v2struct(minPacketDuration, minWAKEEpisodeDuration,...
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

%% Exclude REM from NREM Episodes (ie if the NREM Episode included REM and not just MA)
OrigNREMEpisodes = episodeintervals{2};%this will be a list of acceptable episodes
AllNewNREMEpisodes = [];%will populate this list with snippets of interrupted episodds

REMints = sortrows(REMints);
[~,NREMwRemInside] = InIntervals(REMints(:,1),double(episodeintervals{2}));%NREMs have REM inside
NREMwRemInside(NREMwRemInside == 0) = [];
badNrem = unique(NREMwRemInside);
for nix = length(badNrem):-1:1%backwards bc will be cutting out episodes
    tnrem = badNrem(nix);
    thisNREMint = episodeintervals{2}(tnrem,:);
    OrigNREMEpisodes(tnrem,:) = [];%delete bad episodes from permanent list
    
    whichREMinThisNREM = InIntervals(REMints(:,1),thisNREMint);
    whichREMinThisNREM = find(whichREMinThisNREM);
    for rix = 1:length(whichREMinThisNREM)%for every REM inside this NREM
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

%% Output: buzcode format
% Output 1: 
SleepState.ints.NREMstate = NREMints;
SleepState.ints.REMstate = REMints;
SleepState.ints.WAKEstate = WAKEIntervals;
SleepState.detectorinfo.detectorname = DetectorName;
SleepState.detectorinfo.detectionparams = DetectionParams;
SleepState.detectorinfo.detectiondate = datestr(today,'yyyy-mm-dd');

% Output 2:
SleepStateEpisodes.ints.NREMepisode = episodeintervals{2};
SleepStateEpisodes.ints.REMepisode = episodeintervals{3};
SleepStateEpisodes.ints.WAKEepisode = episodeintervals{1};
SleepStateEpisodes.ints.NREMpacket = packetintervals;
SleepStateEpisodes.ints.MAstate = MAIntervals;
SleepStateEpisodes.ints.MA_REM = MAInREM;
SleepStateEpisodes.detectorinfo.detectionparams.StateDetectionParams = DetectionParams;
SleepStateEpisodes.detectorinfo.detectionparams.EpisodeDetectionParams= detectionparams_episodes;
SleepStateEpisodes.detectorinfo.detector = ['StatesToFinalScoring.m'];
SleepStateEpisodes.detectorinfo.detectiondate = datestr(today,'yyyy-mm-dd');





