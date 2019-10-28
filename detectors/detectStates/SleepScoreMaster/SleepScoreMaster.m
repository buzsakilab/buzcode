function SleepState = SleepScoreMaster(basePath,varargin)
%SleepScoreMaster(basePath,<options>)
%This is the master function for sleep state scoring.
%
%It's strongly recommended that you 
%   1) indicate bad (noisy) channels using bz_getSessionInfo(basePath,'editGUI',true)
%      before running SleepScoreMaster.
%   2) use the 'ignoretime' input to exclude time windows with opto
%       stimulation or behavior with electrical noise
%   3) check the scoring quality using TheStateEditor after running SleepScoreMaster.
%      Use the 'A' key in TheStateEditory to further refine thresholds as
%      needed, and implement sticky thresholds.
%
%INPUT 
%   basePath        folder containing .xml and .lfp files.
%                   basePath and files should be of the form:
%                   'whateverfolder/recordingName/recordingName'
%   (optional)      If no inputs included, select folder(s) containing .lfp
%                   and .xml file in prompt.
%   (optional)      if no .lfp in basePath, option to select multiple 
%                   lfp-containing subfolders
%                          
%   OPTIONS
%   'savedir'       Default: basePath
%   'overwrite'     Overwrite all processing steps (Default: false)
%   'ignoreManual'  Default: false. Overwrite manual scoring from TheStateEditor
%   'savebool'      Default: true. Save anything.
%   'scoretime'     Window of time to score. Default: [0 Inf] 
%                   NOTE: must be continous interval
%   'ignoretime'    Time intervals winthin scoretime to ignore 
%                   (for example, opto stimulation or behavior with artifacts)   
%   'winparms'      [FFT window , smooth window] (Default: [2 15])
%                   (Note: updated from [10 10] based on bimodaility optimization, 6/17/19)
%   'SWWeightsName' Name of file in path (in Dependencies folder) 
%                   containing the weights for the various frequencies to
%                   be used for SWS detection.  
%                   Default is to use Power Spectrum Slope ('PSS'),
%                   If this doesn't work, try 'SWweights.mat' 
%                   or 'SWweightsHPC.mat' for HPC recording.
%   'Notch60Hz'     Boolean 0 or 1.  Value of 1 will notch out the 57.5-62.5 Hz
%                   band, default is 0, no notch.  This can be necessary if
%                   electrical noise.
%   'NotchUnder3Hz' Boolean 0 or 1.  Value of 1 will notch out the 0-3 Hz
%                   band, default is 0, no notch.  This can be necessary
%                   due to poor grounding and low freq movement transients
%   'NotchHVS'      Boolean 0 or 1.  Value of 1 will notch the 4-10 and 
%                   12-18 Hz bands for SW detection, default is 0, no 
%                   notch.  This can be useful in
%                   recordings with prominent high voltage spindles which
%                   have prominent ~16hz harmonics
%   'NotchTheta'    Boolean 0 or 1.  Value of 1 will notch the 4-10 Hz
%                   band for SW detection, default is 0, no notch.  This 
%                   can be useful to
%                   transform the cortical spectrum to approximately
%                   hippocampal, may also be necessary with High Voltage
%                   Spindles
%   'stickytrigger' Implements a "sticky" trigger for SW/EMG threshold 
%                   crossings: metrics must reach halfway between threshold
%                   and opposite peak to count as crossing (reduces
%                   flickering) (default:true)
%   'SWChannels'    A vector list of channels that may be chosen for SW
%                   signal
%   'ThetaChannels' A vector list of channels that may be chosen for Theta
%                   signal
%   'rejectChannels' A vector of channels to exclude from the analysis
%   'saveLFP'       (default:true) to save SleepScoreLFP.lfp.mat file
%   'noPrompts'     (default:false) an option to not prompt user of things
%
%OUTPUT 
%   !THIS IS OUT OF DATE - UPDATE!
%   StateIntervals  structure containing start/end times (seconds) of
%                   NREM, REM, WAKE states and episodes. states is the 
%                   "raw" state scoring. episodes are joined episodes of 
%                   extended (40s) time in a given states, allowing for 
%                   brief interruptions. also contains NREM packets, 
%                   unitary epochs of NREM as described in Watson et al 2016.
%                   saved in a .mat file:
%                   recordingname_SleepScore.mat 
%   
%
% DLevenstein and BWatson 2015/16
%% Recording Selection
%if recname is 'select' or something
%use uigetfile to pick and get list of filenames
%if recname is 'all', get all recordings in a folder and
%then run SleepScoreMaster on each of the filenames'
%if no input arguements... select uigetfile

%Select from no input
if ~exist('basePath','var')
    basePath = uigetdir(cd,...
        'Which recording(s) would you like to state score?');
    if isequal(basePath,0);return;end  
end

%Separate datasetfolder and recordingname
[datasetfolder,recordingname,extension] = fileparts(basePath);
recordingname = [recordingname,extension]; % fileparts parses '.' into extension


%% If there is no .lfp in basePath, choose (multiple?) folders within basePath.
%Select from dataset folder - need to check if .xml/lfp exist
if ~exist(fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file') && ...
    ~exist(fullfile(datasetfolder,recordingname,[recordingname,'.eeg']),'file')
    display(['no .lfp file in basePath, pick a selection of session folders',...
             'containing .lfp files'])
         %Find all basePaths within the topPath
        [basePaths,recordingname] = bz_FindBasePaths(basePath,'select',true);    
end

%If multiple recordings, loop calling SleepScoreMaster with each
numrecs = length(recordingname);
if numrecs > 1 & iscell(recordingname)
    display(['Multiple Recordings (',num2str(numrecs),')'])
    for rr = 1:numrecs
        multibasepath = basePaths{rr};
        SleepScoreMaster(multibasepath,varargin{:})
        close all
    end
    return
elseif numrecs == 1 & iscell(recordingname)
        recordingname = recordingname{1};
end

display(['Scoring Recording: ',recordingname]);

%% inputParse for Optional Inputs and Defaults
p = inputParser;

addParameter(p,'overwrite',false)
addParameter(p,'savebool',true,@islogical)
addParameter(p,'savedir',datasetfolder)
addParameter(p,'scoretime',[0 Inf])
addParameter(p,'ignoretime',[])
addParameter(p,'SWWeightsName','PSS')
addParameter(p,'Notch60Hz',0)
addParameter(p,'NotchUnder3Hz',0)
addParameter(p,'NotchHVS',0)
addParameter(p,'NotchTheta',0)
addParameter(p,'SWChannels',0)
addParameter(p,'ThetaChannels',0)
addParameter(p,'rejectChannels',[]);
addParameter(p,'noPrompts',true);
addParameter(p,'stickytrigger',false);
addParameter(p,'saveLFP',true);
addParameter(p,'winparms',[2 15]);
addParameter(p,'ignoreManual',false)

parse(p,varargin{:})
%Clean up this junk...
overwrite = p.Results.overwrite; 
savedir = p.Results.savedir;
savebool = p.Results.savebool;
scoretime = p.Results.scoretime;
ignoretime = p.Results.ignoretime;
SWWeightsName = p.Results.SWWeightsName;
Notch60Hz = p.Results.Notch60Hz;
NotchUnder3Hz = p.Results.NotchUnder3Hz;
NotchHVS = p.Results.NotchHVS;
NotchTheta = p.Results.NotchTheta;
SWChannels = p.Results.SWChannels;
ThetaChannels = p.Results.ThetaChannels;
rejectChannels = p.Results.rejectChannels;
noPrompts = p.Results.noPrompts;
stickytrigger = p.Results.stickytrigger;
saveLFP = p.Results.saveLFP;
winparms = p.Results.winparms;
ignoreManual = p.Results.ignoreManual; 

%% Database File Management 
savefolder = fullfile(savedir,recordingname);
if ~exist(savefolder,'dir')
    mkdir(savefolder)
end

%Filenames of metadata and SleepState.states.mat file to save
sessionmetadatapath = fullfile(savefolder,[recordingname,'.SessionMetadata.mat']);
%Buzcode outputs
bz_sleepstatepath = fullfile(savefolder,[recordingname,'.SleepState.states.mat']);

%% Check for existing Manual Scoring
if exist(bz_sleepstatepath,'file') && ~overwrite && ~ignoreManual
   SleepState_old = load(bz_sleepstatepath);
   if isfield(SleepState_old.SleepState.detectorinfo,'LastManualUpdate')
       display(['Manual scoring detected... will update the SleepScoreMetrics and AutoScoreInts, ',...
        'but keep previous (manual) scoring.']) 
       display(['To overwrite manual scoring use ''ignoreManual'',true ',...
           'or ''overwrite'',true to overwrite all metrics.'])
       ManScore.ints = SleepState_old.SleepState.ints;
       ManScore.idx = SleepState_old.SleepState.idx;
       ManScore.LastManualUpdate = SleepState_old.SleepState.detectorinfo.LastManualUpdate;
   end
end

%% Get channels not to use
sessionInfo = bz_getSessionInfo(basePath,'noPrompts',noPrompts);
% check that SW/Theta channels exist in rec..
if length(SWChannels) > 1 
    if sum(ismember(SWChannels,sessionInfo.channels)) ~= length(SWChannels)
        error('some of the SW input channels dont exist in this recording...?')
    end   
end
if length(ThetaChannels) > 1 
    if sum(ismember(ThetaChannels,sessionInfo.channels)) ~= length(ThetaChannels)
        error('some of the theta input channels dont exist in this recording...?')
    end   
end

%Is this still needed?/Depreciated?
% if exist(sessionmetadatapath,'file')%bad channels is an ascii/text file where all lines below the last blank line are assumed to each have a single entry of a number of a bad channel (base 0)
%     load(sessionmetadatapath)
%     rejectChannels = [rejectChannels SessionMetadata.ExtracellEphys.BadChannels];
% elseif isfield(sessionInfo,'badchannels')
if isfield(sessionInfo,'badchannels')
    rejectChannels = [rejectChannels sessionInfo.badchannels]; %get badchannels from the .xml
    
end

if isempty(rejectChannels)
    display('No rejected channels - it''s recommended you identify noisy channels to ignore')
end

%% CALCULATE EMG FROM HIGH-FREQUENCY COHERENCE
% Load/Calculate EMG based on cross-shank correlations 
% (high frequency correlation signal = high EMG).  
% Schomburg E.W. Neuron 84, 470?485. 2014)
EMGFromLFP = bz_EMGFromLFP(basePath,'overwrite',overwrite,...
                                     'rejectChannels',rejectChannels,'noPrompts',noPrompts,...
                                     'saveMat',savebool);

%% DETERMINE BEST SLOW WAVE AND THETA CHANNELS
%Determine the best channels for Slow Wave and Theta separation.
%Described in Watson et al 2016, with modifications
SleepScoreLFP = PickSWTHChannel(basePath,...
                            scoretime,SWWeightsName,...
                            Notch60Hz,NotchUnder3Hz,NotchHVS,NotchTheta,...
                            SWChannels,ThetaChannels,rejectChannels,...
                            overwrite,'ignoretime',ignoretime,...
                            'noPrompts',noPrompts,'saveFiles',saveLFP&savebool,...
                            'window',winparms(1),'smoothfact',winparms(2),'IRASA',true);

%% CLUSTER STATES BASED ON SLOW WAVE, THETA, EMG

%Calculate the scoring metrics: broadbandLFP, theta, EMG
display('Quantifying metrics for state scoring')
[SleepScoreMetrics,StatePlotMaterials] = ClusterStates_GetMetrics(...
                                           basePath,SleepScoreLFP,EMGFromLFP,overwrite,...
                                           'onSticky',stickytrigger,'ignoretime',ignoretime,...
                                           'window',winparms(1),'smoothfact',winparms(2),...
                                           'IRASA',true,'ThIRASA',true);
                                       
%Use the calculated scoring metrics to divide time into states
display('Clustering States Based on EMG, SW, and TH LFP channels')
[ints,idx,MinTimeWindowParms] = ClusterStates_DetermineStates(SleepScoreMetrics);


                                
%% RECORD PARAMETERS from scoring
detectionparms.userinputs = p.Results;
detectionparms.MinTimeWindowParms = MinTimeWindowParms;
detectionparms.SleepScoreMetrics = SleepScoreMetrics;

% note and keep special version of original hists and threshs

SleepState.ints = ints;
SleepState.idx = idx;
SleepState.detectorinfo.detectorname = 'SleepScoreMaster';
SleepState.detectorinfo.detectionparms = detectionparms;
SleepState.detectorinfo.detectionparms.histsandthreshs_orig = detectionparms.SleepScoreMetrics.histsandthreshs;
SleepState.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');
SleepState.detectorinfo.StatePlotMaterials = StatePlotMaterials;

%Put old manual scoring back in
if exist('ManScore','var')
    SleepState.AutoScoreInts = SleepState.ints;
    SleepState.ints = ManScore.ints;
    SleepState.idx = ManScore.idx;
    SleepState.detectorinfo.LastManualUpdate = ManScore.LastManualUpdate;
end

%Saving SleepStates
save(bz_sleepstatepath,'SleepState');

%% MAKE THE STATE SCORE OUTPUT FIGURE
%ClusterStates_MakeFigure(stateintervals,stateIDX,figloc,SleepScoreMetrics,StatePlotMaterials);
try
    ClusterStates_MakeFigure(SleepState,basePath,noPrompts);
    disp('Figures Saved to StateScoreFigures')
catch
    disp('Figure making error')
end

%% JOIN STATES INTO EPISODES

% Extract states, Episodes, properly organize params etc, prep for final saving
display('Calculating/Saving Episodes')
StatesToEpisodes(SleepState,basePath);

display(['Sleep Score ',recordingname,': Complete!']);

%% PROMPT USER TO MANUALLY CHECK DETECTION WITH THESTATEEDITOR
if ~noPrompts
    str = input('Would you like to check detection with TheStateEditor? [Y/N] ','s');
    switch str
        case {'Y','y',''}
            TheStateEditor([basePath,filesep,recordingname])
        case {'N','n'}
        otherwise
            display('Unknown input..... you''ll have to load TheStateEditor on your own')
    end
end

end

