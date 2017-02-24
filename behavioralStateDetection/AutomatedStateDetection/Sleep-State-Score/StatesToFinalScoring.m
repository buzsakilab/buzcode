function StateIntervals = StatesToFinalScoring(NREMints,WAKEints,REMints)

minPacketDuration = 30;
maxMicroarousalDuration = 100;

maxEpisodeDuration = 40;
minSWSEpisodeDuration = 20;
minWAKEEpisodeDuration = 20;
minREMEpisodeDuration = 20;


SWSlengths = NREMints(:,2)-NREMints(:,1);
packetintervals = NREMints(SWSlengths>=minPacketDuration,:);

WAKElengths = WAKEints(:,2)-WAKEints(:,1);
MAIntervals = WAKEints(WAKElengths<=maxMicroarousalDuration,:);
WAKEIntervals = WAKEints(WAKElengths>maxMicroarousalDuration,:);

[episodeintervals{2}] = IDStateEpisode(NREMints,maxEpisodeDuration,minSWSEpisodeDuration);
[episodeintervals{1}] = IDStateEpisode(WAKEints,maxEpisodeDuration,minWAKEEpisodeDuration);
[episodeintervals{3}] = IDStateEpisode(REMints,maxEpisodeDuration,minREMEpisodeDuration);


%% Save
StateIntervals.NREMstate = NREMints;
StateIntervals.REMstate = REMints;
StateIntervals.WAKEstate = WAKEIntervals;
StateIntervals.NREMepisode = episodeintervals{2};
StateIntervals.REMepisode = episodeintervals{3};
StateIntervals.WAKEeposode = episodeintervals{1};
StateIntervals.NREMpacket = packetintervals;
StateIntervals.MAstate = MAIntervals;



%% Identify MAs preceeding REM
% [ ~,~,MA_REM,~ ] = FindIntsNextToInts(StateIntervals.MAstate,StateIntervals.REMstate);
% realMA = setdiff(1:length(StateIntervals.MAstate(:,1)),MA_REM);
% StateIntervals.MA_REM = StateIntervals.MAstate(MA_REM,:);
% StateIntervals.MAstate = StateIntervals.MAstate(realMA,:);

