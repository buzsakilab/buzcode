function [SleepStateEpisodes] = StatesToEpisodes(SleepState,basePath)
% Takes a buzcode structure of sleep state data and imposes max and
% min durations and interruption criteria to make WAKE, NREM, REM Episodes,
% Packets, Microarousals and Microarusals in REM, as defined in 
% Watson et al 2017. Returns the SleepStateEpisodes in buzcode format
% along with any detectionparms from initial state scoring. Can save as
% states.mat file in basePath
% structures.
%
% INPUTS
%   SleepState  buzcode states.mat structure. 
%               Expected to have fields:
%                   SleepState.ints.NREMstate
%                   SleepState.ints.WAKEstate
%                   SleepState.ints.REMstate
%               Hopefully has fields:
%                   SleepState.detectorinfo
%   basePath    (optional) a basepath to save the output
% 
% OUTPUTS
%	SleepStateEpisodes - Struct array of relatively more processed epochs
%   with minimum durations and maximum interruptions imposed.  Fields below
%       .ints.WAKEpisodes  WAKEints inputs with at least [minWAKEEpisodeDuration] 
%                          long and without interruptions more than [maxWAKEEpisodeInterruption]
%       .ints.NREMepisodes
%       .ints.REMepisodes
%       .ints.NREMpackets
%       .ints.MAstate
%       .ints.MAInREM 
%       .detectorinfo      a structure with detector information
% 

%
% Called in [SleepScoreMaster] and [TheStateEditor "saveStates" & "ReClusterStates_In" functions]  
% Dan Levenstein & Brendon Watson 2016
%%

if ~isfield(SleepState.ints,'NREMstate')
    SleepState.ints.NREMstate = zeros(0,2);
end
if ~isfield(SleepState.ints,'WAKEstate')
    SleepState.ints.WAKEstate = zeros(0,2);
end
if ~isfield(SleepState.ints,'REMstate')
    SleepState.ints.REMstate = zeros(0,2);
end

%Pull inputs from SleepState
NREMints =  SleepState.ints.NREMstate;
WAKEints =  SleepState.ints.WAKEstate;
REMints =  SleepState.ints.REMstate;

%Duration parameters
minPacketDuration = 30;
minWAKEEpisodeDuration = 20;
minNREMEpisodeDuration = 20;
minREMEpisodeDuration = 20;

maxMicroarousalDuration = 100;

maxWAKEEpisodeInterruption = 40;
maxNREMEpisodeInterruption = maxMicroarousalDuration;
maxREMEpisodeInterruption = 40;

detectionparms_episodes = v2struct(minPacketDuration, minWAKEEpisodeDuration,...
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
% Output: SleepStateEpisodes
SleepStateEpisodes.ints.NREMepisode = episodeintervals{2};
SleepStateEpisodes.ints.REMepisode = episodeintervals{3};
SleepStateEpisodes.ints.WAKEepisode = episodeintervals{1};
SleepStateEpisodes.ints.NREMpacket = packetintervals;
SleepStateEpisodes.ints.MA = MAIntervals;
SleepStateEpisodes.ints.MA_REM = MAInREM;
SleepStateEpisodes.detectorinfo.originaldetectorinfo = SleepState.detectorinfo;
SleepStateEpisodes.detectorinfo.detectionparms.EpisodeDetectionParms= detectionparms_episodes;
SleepStateEpisodes.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');


%% Save

%Saving SleepStateEpisodes (most interpreted/processed) - bzStyle
if exist('basePath','var')
    baseName = bz_BasenameFromBasepath(basePath);
    bz_sleepstateepisodespath = fullfile(basePath,[baseName,'.SleepStateEpisodes.states.mat']);
    save(bz_sleepstateepisodespath,'SleepStateEpisodes');
end



