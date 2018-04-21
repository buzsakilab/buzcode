function [INT, IDX, t_IDX] = ClusterStates_DetermineStates(SleepScoreMetrics,MinWinParams,histsandthreshs)
% can input histsandthreshs from externally if needed... ie via manual
% selection in stateeditor

%% Basic parameters
% Min Win Parameters (s)
if exist('MinWinParams','var')
     v2struct(MinWinParams)
% 	fn = fieldnames(MinWinParams);
%     for a = 1:length(fn);
%         eval([fn{a} '=MinWinParams.' fn{a} ';']);
%     end
else%defaults as follows:
    minSWS = 6;
    minWnexttoREM = 6;
    minWinREM = 6;       
    minREMinW = 6;
    minREM = 6;
    minWAKE = 6;
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
NREMtimes = (broadbandSlowWave >swthresh);
MOVtimes = (broadbandSlowWave<swthresh & EMG>EMGthresh);
if THthresh ~= 0
    REMtimes = (broadbandSlowWave<swthresh & EMG<EMGthresh & thratio>THthresh);
else % THthresh = 0;
    REMtimes =(broadbandSlowWave<swthresh & EMG<EMGthresh);
end    

%%
%OLD:
%Index Vector: SWS=2, REM=3, MOV=6, NonMOV=1.   
%(Separate MOV for REM, then join later)
%IDX = SWStimes+2*REMtimes+5*MOVtimes+1;

%NEW: No separation of MOV and NonMOV WAKE
%Index Vector: NREM=2, REM=3, WAKE=1. 
IDX = NREMtimes+2*REMtimes+1;

%Start/end offset due to FFT


%% Minimum Interuptions
INT = IDXtoINT_ss(IDX,3);


%Make the following repeated chunks of code into a single function.

%SWS  (to NonMOV)
Sints = INT{2};
Slengths = Sints(:,2)-Sints(:,1);
shortSints = {Sints(find(Slengths<=minSWS),:)};
shortSidx = INTtoIDX_ss(shortSints,length(IDX));
%Change Short SWS to Wake
IDX(shortSidx==1) = 1;   
INT = IDXtoINT_ss(IDX,3);

%NonMOV next to REM   (to REM)
Wints = INT{1};
trans = (diff(IDX)); %All State Transitions
WRtrans = find((trans)==-2)+1;  %Just transitions between WAKE, REM
%Convert to interval indices
[~,WRtransON] = intersect(Wints(:,1),WRtrans);
WRtrans = find((trans)==2);
[~,WRtransOFF] = intersect(Wints(:,2),WRtrans);
WRtrans = union(WRtransON,WRtransOFF); %On or offset are RW
%Find WAKE intervals that border REM and are less than min
Wlengths = Wints(:,2)-Wints(:,1);
shortWRints = find(Wlengths(WRtrans)<=minWnexttoREM);
shortWRints = WRtrans(shortWRints);
shortWRints = {Wints(shortWRints,:)};
shortWRidx = INTtoIDX_ss(shortWRints,length(IDX));
%Convert wake to rem
IDX(shortWRidx==1) = 3;
INT = IDXtoINT_ss(IDX,3);


%NonMOV in REM   (to REM)
Wints = INT{1};
trans = (diff(IDX)); %All State Transitions
WRtrans = find((trans)==-2)+1;  %Just transitions between WAKE, REM
%Convert to interval indices
[~,WRtransON] = intersect(Wints(:,1),WRtrans);
WRtrans = find((trans)==2);
[~,WRtransOFF] = intersect(Wints(:,2),WRtrans);
WRtrans = intersect(WRtransON,WRtransOFF); %Both onset and offset are RW
%Find WAKE intervals that border REM and are less than min
Wlengths = Wints(:,2)-Wints(:,1);
shortWRints = find(Wlengths(WRtrans)<=minWinREM);
shortWRints = WRtrans(shortWRints);
shortWRints = {Wints(shortWRints,:)};
shortWRidx = INTtoIDX_ss(shortWRints,length(IDX));
%Convert wake to rem
IDX(shortWRidx==1) = 3;
IDX(IDX==6) = 1; %Convert NonMOV to WAKE
INT = IDXtoINT_ss(IDX,3);


%REM in WAKE   (to WAKE)
Rints = INT{3};
trans = (diff(IDX)); %All State Transitions
WRtrans = find((trans)==2)+1;  %Just transitions between WAKE, REM
%Convert to interval indices
[~,WRtransON] = intersect(Rints(:,1),WRtrans);
WRtrans = find((trans)==-2);
[~,WRtransOFF] = intersect(Rints(:,2),WRtrans);
WRtrans = intersect(WRtransON,WRtransOFF); %Both onset and offset are RW
%Find WAKE intervals that border REM and are less than min
Rlengths = Rints(:,2)-Rints(:,1);
shortWRints = find(Rlengths(WRtrans)<=minREMinW);
shortWRints = WRtrans(shortWRints);
shortWRints = {Rints(shortWRints,:)};
shortWRidx = INTtoIDX_ss(shortWRints,length(IDX));
%Convert REM to WAKE
IDX(shortWRidx==1) = 1;
INT = IDXtoINT_ss(IDX,3);

%REM (only applies to REM in the middle of SWS)    (to WAKE)
Rints = INT{3};
Rlengths = Rints(:,2)-Rints(:,1);
shortRints = {Rints(find(Rlengths<=minREM),:)};
shortRidx = INTtoIDX_ss(shortRints,length(IDX));

IDX(shortRidx==1) = 1;
INT = IDXtoINT_ss(IDX,3);


%WAKE   (to SWS)     essentiall a minimum MA time
Wints = INT{1};
Wlengths = Wints(:,2)-Wints(:,1);
shortWints = {Wints(find(Wlengths<=minWAKE),:)};
shortWidx = INTtoIDX_ss(shortWints,length(IDX));
IDX(shortWidx==1) = 2;

INT = IDXtoINT_ss(IDX,3);

%SWS  (to NonMOV)
Sints = INT{2};
Slengths = Sints(:,2)-Sints(:,1);
shortSints = {Sints(find(Slengths<=minSWS),:)};
shortSidx = INTtoIDX_ss(shortSints,length(IDX));
%Change Short SWS to Wake
IDX(shortSidx==1) = 1;   
INT = IDXtoINT_ss(IDX,3);




%% Pad time to match recording time
offset = SleepScoreMetrics.t_clus(1)-1; %t_FFT(1)-1;

INT = cellfun(@(x) x+offset,INT,'UniformOutput',false);

IDX = INTtoIDX_ss(INT,reclength);
t_IDX = 1:length(IDX);
IDX = IDX';

end

