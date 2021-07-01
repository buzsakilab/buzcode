function [ints, idx, MinTimeWindowParms] = ClusterStates_DetermineStates(...
    SleepScoreMetrics,MinTimeWindowParms,histsandthreshs)
% can input histsandthreshs from externally if needed... ie via manual
% selection in stateeditor

%% Basic parameters
% Min Win Parameters (s)
if exist('MinTimeWindowParms','var') && ~isempty(MinTimeWindowParms)
     v2struct(MinTimeWindowParms)
else%defaults as follows:
    minSWSsecs = 6;
    minWnexttoREMsecs = 6;
    minWinREMsecs = 6;       
    minREMinWsecs = 6;
    minREMsecs = 6;
    minWAKEsecs = 6;
    MinTimeWindowParms = v2struct(minSWSsecs,minWnexttoREMsecs,minWinREMsecs,...
        minREMinWsecs,minREMsecs,minWAKEsecs);
end

% handling variables for determining thresholds/cutoffs
if exist('histsandthreshs','var')
    hat2 = histsandthreshs;%bc will overwrite below
    v2struct(SleepScoreMetrics)
    histsandthreshs = hat2;
else
    v2struct(SleepScoreMetrics)
end
v2struct(histsandthreshs)%Expand and get values out of these fields

%% Re-Do this code (should be same as in ClusterStates_GetParams.m) to see if theta is bimodal

%This switch turns on a "schmidt trigger", or sticky trigger,
%which means that threshold crossings have to reach the
%midpoint between the dip and the opposite peak, this
%reduces noise. Passed through via histsandthreshs from checkboxes in
%TheStateEditor or 'stickytrigger',true in SleepScoreMaster via GetMetrics
if ~exist('stickySW','var'); stickySW = false; end
if ~exist('stickyTH','var'); stickyTH = false; end
if ~exist('stickyEMG','var'); stickyEMG = false; end


[~,~,~,~,NREMtimes] = bz_BimodalThresh(broadbandSlowWave(:),'startbins',15,...
    'setthresh',swthresh,'diptest',false,'Schmidt',stickySW,'0Inf',true);

[~,~,~,~,hightheta] = bz_BimodalThresh(thratio(:),'startbins',15,...
    'setthresh',THthresh,'diptest',false,'Schmidt',stickyTH,'0Inf',true);

[~,~,~,~,highEMG] = bz_BimodalThresh(EMG(:),'startbins',15,...
    'setthresh',EMGthresh,'diptest',false,'Schmidt',stickyEMG,'0Inf',true);

REMtimes = (~NREMtimes & ~highEMG & hightheta);

%ACTIVE/QUIET WAKE:
WAKEtimes = ~NREMtimes & ~REMtimes;
QWAKEtimes =  WAKEtimes & ~hightheta; %Used later if QWake scored

%%
%Start/end offset due to FFT

%Construct IDX STRUCTURE FOR bz_IDXtoINT
IDX.statenames = {'WAKE','','NREM','','REM'};
IDX.timestamps = t_clus; %Timestamps pulled from clustering (in ClusterStates_GetMetrics)
IDX.states = zeros(size(IDX.timestamps));
IDX.states(NREMtimes) = 3;
IDX.states(REMtimes) = 5;
IDX.states(WAKEtimes) = 1;

if ~exist('scoreQW','var'); scoreQW = false; end
if scoreQW
    IDX.states(QWAKEtimes) = 2;
    IDX.statenames{2} = 'QWAKE';
end
%% Fill in gaps with 0, by converting to intervals and then back to indices/timestamps
INT = bz_IDXtoINT(IDX);
IDX = bz_INTtoIDX(INT,'statenames',{'WAKE','','NREM','','REM'});
%% Minimum Interuptions



%Make the following repeated chunks of code into a single function.

%Short NREM -> WAKE
Sdur = diff(INT.NREMstate,[],2);
shortSints = Sdur<=minSWSsecs;
shortSidx = InIntervals(IDX.timestamps,INT.NREMstate(shortSints,:));
IDX.states(shortSidx) = 1;   
INT = bz_IDXtoINT(IDX);

%Short WAKE (next to REM) -> REM
Wdur = diff(INT.WAKEstate,[],2);
shortWints = Wdur<=minWnexttoREMsecs;
[~,~,WRints,~] = FindIntsNextToInts(INT.WAKEstate,INT.REMstate,1);
[~,~,~,RWints] = FindIntsNextToInts(INT.REMstate,INT.WAKEstate,1);
shortWRidx = InIntervals(IDX.timestamps,INT.WAKEstate(shortWints & (WRints | RWints),:));
IDX.states(shortWRidx) = 5;   
INT = bz_IDXtoINT(IDX);

%Short REM (in WAKE) -> WAKE
Rdur = diff(INT.REMstate,[],2);
shortRints =Rdur<=minREMinWsecs;
[~,~,~,WRints] = FindIntsNextToInts(INT.WAKEstate,INT.REMstate,1);
[~,~,RWints,~] = FindIntsNextToInts(INT.REMstate,INT.WAKEstate,1);
shortRWidx = InIntervals(IDX.timestamps,INT.REMstate(shortRints & (WRints & RWints),:));
IDX.states(shortRWidx) = 1;   
INT = bz_IDXtoINT(IDX);

%Remaining Short REM (in NREM) -> WAKE
Rdur = diff(INT.REMstate,[],2);
shortRints = Rdur<=minREMsecs;
shortRidx = InIntervals(IDX.timestamps,INT.REMstate(shortRints,:));
IDX.states(shortRidx) = 1;   
INT = bz_IDXtoINT(IDX);

%WAKE   (to SWS)     essentiall a minimum MA time
Wdur = diff(INT.WAKEstate,[],2);
shortWints = Wdur<=minWAKEsecs;
shortWidx = InIntervals(IDX.timestamps,INT.WAKEstate(shortWints,:));
IDX.states(shortWidx) = 3;   
INT = bz_IDXtoINT(IDX);

%SWS  (to NonMOV)
Sdur = diff(INT.NREMstate,[],2);
shortSints = Sdur<=minSWSsecs;
shortSidx = InIntervals(IDX.timestamps,INT.NREMstate(shortSints,:));
IDX.states(shortSidx) = 1;   
INT = bz_IDXtoINT(IDX);

%% Output... lol

ints = INT;
idx = IDX;



end

