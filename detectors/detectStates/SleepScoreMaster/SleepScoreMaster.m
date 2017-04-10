function SleepScoreMaster(basePath,varargin)
%SleepScoreMaster(datasetfolder,recordingname)
%This is the master function for sleep state scoring.
%
%INPUT 
%   basePath        folder containing .xml and .lfp files.
%                   basePath and files should be of the form:
%                   whateverfolder/recordingName/recordingName.lfp
%   (optional)      If no inputs included, select folder(s) containing .lfp
%                   and .xml file in prompt.
%   (optional)      if no .lfp in basePath, option to select multiple 
%                   lfp-containing subfolders
%                          
%   OPTIONS
%   'savedir'       Default: datasetfolder
%   'overwrite'     Default: false
%   'savebool'      Default: true
%   'scoretime'     Default: [0 Inf]
%   'badchannels'   file datasetfolder/recordingname/'bad_channels.txt'
%                   that lists channels will omit certain channels from EMG
%                   detection and LFP selection
%   'SWWeightsName' Name of file in path (in Dependencies folder) 
%                   containing the weights for the various frequencies to
%                   be used for SWS detection.  Default is 'SWweights.mat'
%                     - For hippocampus-only recordings, enter
%                     'SWWeightsHPC.mat' for this
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
%   'SWChannels'    A vector list of channels that may be chosen for SW
%                   signal
%   'ThetaChannels' A vector list of channels that may be chosen for Theta
%                   signal
%
%OUTPUT
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

%% Parameter setting
% Min Win Parameters (s): basic detection paramaters (seconds)
minSWS = 6;
minWnexttoREM = 6;
minWinREM = 6;       
minREMinW = 6;
minREM = 6;
minWAKE = 6;
MinWinParams = v2struct(minSWS,minWnexttoREM,minWinREM,minREMinW,minREM,minWAKE);

%% Recording Selection
%if recname is 'select' or something
%use uigetfile to pick and get list of filenames
%if recname is 'all', get all recordings in a folder and
%then run SleepScoreMaster on each of the filenames'
%if no input arguements... select uigetfile

%Select from no input
if ~exist('basePath','var')
    DIRECTORYNAME = uigetdir(cd,...
        'Which recording(s) would you like to state score?');
    if isequal(DIRECTORYNAME,0);return;end  
    [datasetfolder,recordingname] = fileparts(DIRECTORYNAME); 
else
    %Separate datasetfolder and recordingname
    [datasetfolder,recordingname] = fileparts(basePath);
end


if ~exist('SWWeightsName','var')
    SWWeightsName = 'SWweights.mat';
end



%% If there is no .lfp in basePath, choose (multiple?) folders within basePath.
%Select from dataset folder - need to check if .xml/lfp exist
if ~exist(fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file')
        foldercontents = dir(basePath);
        possiblerecordingnames = {foldercontents([foldercontents.isdir]==1).name};
        [s,v] = listdlg('PromptString','Which recording(s) would you like to state score?',...
                        'ListString',possiblerecordingnames);
        recordingname = possiblerecordingnames(s);
end

%If multiple recordings, loop
numrecs = length(recordingname);
if numrecs > 1 & iscell(recordingname)
    display(['Multiple Recordings (',num2str(numrecs),')'])
    for rr = 1:numrecs
        multibasepath = fullfile(basePath,recordingname{rr});
        SleepScoreMaster(multibasepath,varargin{:})
        close all
    end
    return
elseif numrecs == 1 & iscell(recordingname)
        recordingname = recordingname{1};
end

display(['Scoring Recording: ',recordingname]);

%% Deal with input options from varargin
%none yet, but will do this with inputParser when we have input options

%possible variable input options
%'timewin' - only state score a subset of the recording
%'HPCsites' - site indices for HPC probes - will only check these for theta
%           if applicable
%'figloc' - secondardy folder to save figures to
%'spikegroups' - if not in the .xml file
%'SWChannel', 'ThetaChannel' - can enter manually instead of determining
%                               algorithmically

%% inputParse for Optional Inputs and Defaults
p = inputParser;

defaultOverwrite = false;    %Pick new and Overwrite existing ThLFP, SWLFP?
defaultSavebool = true;    %Save Stuff (EMG, LFP)

defaultSavedir = datasetfolder;

defaultScoretime = [0 Inf];
defaultSWWeightsName = 'SWweights.mat';
defaultNotch60Hz = 0;
defaultNotchUnder3Hz = 0;
defaultNotchHVS = 0;
defaultNotchTheta = 0;
defaultSWChannels = 0;
defaultThetaChannels = 0;

addParameter(p,'overwrite',defaultOverwrite)
addParameter(p,'savebool',defaultSavebool,@islogical)
addParameter(p,'savedir',defaultSavedir)
addParameter(p,'scoretime',defaultScoretime)
addParameter(p,'SWWeightsName',defaultSWWeightsName)
addParameter(p,'Notch60Hz',defaultNotch60Hz)
addParameter(p,'NotchUnder3Hz',defaultNotchUnder3Hz)
addParameter(p,'NotchHVS',defaultNotchHVS)
addParameter(p,'NotchTheta',defaultNotchTheta)
addParameter(p,'SWChannels',defaultNotchTheta)
addParameter(p,'ThetaChannels',defaultNotchTheta)

parse(p,varargin{:})
%Clean up this junk...
overwrite = p.Results.overwrite; 
savebool = p.Results.savebool;
savedir = p.Results.savedir;
scoretime = p.Results.scoretime;
SWWeightsName = p.Results.SWWeightsName;
Notch60Hz = p.Results.Notch60Hz;
NotchUnder3Hz = p.Results.NotchUnder3Hz;
NotchHVS = p.Results.NotchHVS;
NotchTheta = p.Results.NotchTheta;
SWChannels = p.Results.SWChannels;
ThetaChannels = p.Results.ThetaChannels;

%% Database File Management 
savefolder = fullfile(savedir,recordingname);

if ~exist(savefolder,'dir')
    mkdir(savefolder)
end

%Figure locations
figloc = [fullfile(savefolder,'StateScoreFigures'),'/'];
if ~exist(figloc,'dir')
    mkdir(figloc)
end


%Filenames for EMG, thLFP, and swLFP .mat files in the database.
%ALL OF THESE NEED TO BE CLEANED INTO BUZCODE FORMAT
%Theta/SWLFP are depreciated - replaced with scorelfppath
sessionmetadatapath = fullfile(savefolder,[recordingname,'_SessionMetadata.mat']);
scorelfppath = fullfile(savefolder,[recordingname,'_SleepScoreLFP.mat']);
EMGpath = fullfile(savefolder,[recordingname '_EMGCorr.mat']);
%Filenames for State and Event .mat files.
sleepstatepath = fullfile(savefolder,[recordingname,'_SleepScore.mat']);
%Filenames for StateCluster Metrics (broadband/theta)
scoremetricspath = fullfile(savefolder,[recordingname,'_StateScoreMetrics.mat']);

%Buzcode outputs
bz_sleepstatepath = fullfile(savefolder,[recordingname,'.SleepState.states.mat']);
bz_scorelfppath = fullfile(savefolder,[recordingname,'.SleepScoreLFP.LFP.mat']);
bz_EMGpath = fullfile(savefolder,[recordingname '.EMGCorr.LFP.mat']);
bz_scoremetricspath = fullfile(savefolder,[recordingname,'.SleepScoreMetrics.mat']);


%Filename for .lfp file
if exist (fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.lfp']);
elseif exist (fullfile(datasetfolder,recordingname,[recordingname,'.lfp']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.lfp']);
elseif exist (fullfile(datasetfolder,recordingname,[recordingname,'.eeg']),'file')
    rawlfppath = fullfile(datasetfolder,recordingname,[recordingname,'.eeg']);
    
elseif ~overwrite
    display('No .lfp file... but using saved files so maybe it''s ok!')
else
    display('No .lfp file')
    return
end

%% Get channels not to use
if exist(sessionmetadatapath,'file')%bad channels is an ascii/text file where all lines below the last blank line are assumed to each have a single entry of a number of a bad channel (base 0)
    load(sessionmetadatapath)
    rejectchannels = SessionMetadata.ExtracellEphys.BadChannels;
else
    rejectchannels = [];
end


%% CALCULATE EMG FROM HIGH-FREQUENCY COHERENCE
% Load/Calculate EMG based on cross-shank correlations 
% (high frequency correlation signal = high EMG).  
% Schomburg E.W. Neuron 84, 470?485. 2014)
% Do this before lfp load because finding proper lfps will depend on this.
% If EMG is already calculated and in it's own .mat, then load, otherwise
% calculate this

% sf_EMG = 2;
if ~exist(EMGpath,'file') || overwrite;
%     [PATHSTR] = fileparts([datasetfolder '/' recordingname]);
%     if exist(fullfile(PATHSTR,'bad_channels.txt'),'file')%bad channels is an ascii/text file where all lines below the last blank line are assumed to each have a single entry of a number of a bad channel (base 0)
%         t = ReadBadChannels_ss(PATHSTR);
%         rejectchannels = cat(1,rejectchannels(:),t(:));
%     end % this should be replaced by a search of the meta data file instead of a separate bad_chans file

    display('Calculating EMG')
    %Run the EMG Extraction 
    bz_EMGFromLFP(rawlfppath,'restrict',scoretime,[]...
                                     ,'rejectChannels',rejectchannels);
%     if savebool
%         %Old Format - update to buzcode
%         save(EMGpath,'EMGCorr','sf_EMG')
%         
%         %Buzcode format - wrap this into EMGCorrForSleepScore - make a
%         %standalone buzcode detector
%         CorrEMG.data = EMGCorr;
%         CorrEMG.sf = sf_EMG;
%         save(bz_EMGpath,'CorrEMG')
%     end

else
    display('EMG aleady calculated: Loading...')

end

    load(EMGpath,'EMGCorr')
% EMG = EMGCorr(:,2);
% clear EMGCorr


%% DETERMINE BEST SLOW WAVE AND THETA CHANNELS
if ~exist(scorelfppath,'file') || overwrite; % if no lfp file already, load lfp and make lfp file?

    display('Picking SW and TH Channels')
    [SleepScoreLFP] = PickSWTHChannel(datasetfolder,recordingname,figloc,scoretime,SWWeightsName,Notch60Hz,NotchUnder3Hz,NotchHVS,NotchTheta,SWChannels,ThetaChannels,rejectchannels);
    if savebool
        save(bz_scorelfppath,'SleepScoreLFP');  
    end
else
    display('SlowWave and Theta LFP Channels Already Extracted, Loading...')
    load(bz_scorelfppath)
        
end

% %CAN THIS BE REMOVED?
% if ~exist('SWfreqlist','var')
%         load(SWWeightsName)%load default weights which would have been used for these older scorings... so they can be saved
% end


%% CLUSTER STATES BASED ON SLOW WAVE, THETA, EMG

display('Quantifying metrics for state scoring')
[broadbandSlowWave,thratio,EMG_env,t_EMG,t_clus,badtimes,reclength,histsandthreshs,...
    swFFTfreqs,swFFTspec,thFFTfreqs,thFFTspec] = ClusterStates_GetMetrics(SleepScoreLFP,EMGCorr);

display('Clustering States Based on EMG, SW, and TH LFP channels')
[stateintervals,stateIDX,~] = ClusterStates_DetermineStates(...
    broadbandSlowWave,thratio,t_clus,EMG_env,histsandthreshs,MinWinParams,reclength,...
    figloc);

ClusterStates_MakeFigure(stateintervals,stateIDX,figloc,swFFTfreqs,swFFTspec,...
    thFFTfreqs,thFFTspec,t_clus,recordingname,broadbandSlowWave,thratio,EMG_env,t_EMG);

if savebool
    %Should save (downsampled to what's used in clusterstates...)
    %sw/thLFP in scoremetricspath here!
    save(scoremetricspath,...
        'broadbandSlowWave','thratio','EMG_env','t_clus',...
        'SWchannum','THchannum','badtimes','reclength','histsandthreshs',...
        'SWfreqlist','SWweights','SWWeightsName','Notch60Hz',...
        'NotchUnder3Hz','NotchHVS','NotchTheta')
    
    save(bz_scoremetricspath,'SleepScoreMetrics')
    
end

%% JOIN STATES INTO EPISODES

NREMints = stateintervals{2};
REMints = stateintervals{3};
WAKEints = stateintervals{1};

[StateIntervals,SleepState,durationprams] = StatesToFinalScoring(NREMints,WAKEints,REMints);

%Old Style - remove once bzCompadible with StateEditor
StateIntervals.metadata.SWchannum = SWchannum;
StateIntervals.metadata.THchannum = THchannum;
save(sleepstatepath,'StateIntervals');

%bzStyle
SleepState.detectorparms.SWchannum = SWchannum;
SleepState.detectorparms.THchannum = THchannum;
SleepState.detectorparms.durationprams = durationprams;
SleepState.detectorname = 'SleepScoreMaster';
SleepState.detectiondate = today;
save(bz_sleepstatepath,'SleepState');

display(['Sleep Score ',recordingname,': Complete!']);

end

