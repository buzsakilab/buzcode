function ConvertTwoLFPFilessToOne(datasetfolder,recordingname,varargin)
% Converts previous _SWLFP.mat and ThetaLFP.mat files to a single
% StateScoreLFP.mat file.  
%INPUT (optional)   If no inputs included, select folder containing .lfp
%                   and .xml file in prompt.
%                   
%   datasetfolder   Top level folder in which the dataset resides. 
%                   For example:
%                   '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/'
%                   -if not included, prompt comes up to
%                   navigate to the folder holding your recording.
%   recordingname   (optional)
%                   Name of the recording, this will be the name of the
%                   folder in which the .lfp file and other files reside.
%                   For example, the .lfp file should be:
%                   'datasetfolder/recordingname/recordingname.lfp'
%                   ... it is also assumed that this serves as the basename
%                   for the files for instance data will be at
%                   /datasetfolder/recordingname/recordingname.lfp
%   'savedir'       Default: datasetfolder
%
% Otherinputs not used
%
% Brendon Watson 2016

if ~exist('datasetfolder','var')
    DIRECTORYNAME = uigetdir('',...
        'Which recording(s) would you like to state score?');
    if isequal(DIRECTORYNAME,0);return;end  
    [datasetfolder,recordingname] = fileparts(DIRECTORYNAME); 
end
  
%Select from dataset folder
switch recordingname
    case 'select'
        foldercontents = dir(datasetfolder);
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
        SleepScoreMaster(datasetfolder,recordingname{rr},varargin{:})
        close all
    end
    return
elseif numrecs == 1 & iscell(recordingname)
        recordingname = recordingname{1};
end

display(['Scoring Recording: ',recordingname]);

%% inputParse for Optional Inputs and Defaults
p = inputParser;

defaultOverwrite = false;    %Pick new and Overwrite existing ThLFP, SWLFP?
defaultSavebool = true;    %Save Stuff (EMG, LFP)
defaultSpindledelta = false; %Detect spindles/delta?

defaultSavedir = datasetfolder;

defaultScoretime = [0 Inf];

addParameter(p,'overwrite',defaultOverwrite,@islogical)
addParameter(p,'savedir',defaultSavedir)


parse(p,varargin{:})
%Clean up this junk...
savedir = p.Results.savedir;

%% Database File Management 
savefolder = fullfile(savedir,recordingname);
scorelfppath = fullfile(savefolder,[recordingname,'_SleepScoreLFP.mat']);
thetalfppath = fullfile(savefolder,[recordingname,'_ThetaLFP.mat']);
swlfppath = fullfile(savefolder,[recordingname,'_SWLFP.mat']);

if ~exist(savefolder,'dir')
    mkdir(savefolder)
end

%% Convert old double-files to new single file
if ((~exist(thetalfppath,'file') && ~exist(swlfppath,'file')) && ~exist(scorelfppath,'file'))  % if no lfp file already, load lfp and make lfp file?
    display('No SW and TH Channels found, skipping')
else
    display('SW and TH Channels Already Extracted, Loading...')
    %For updating state score LFP storage...
    if ~exist(scorelfppath,'file')
        load(swlfppath,'swLFP','SWchannum','sf_LFP')
        load(thetalfppath,'thLFP','THchannum','sf_LFP')
        if sf_LFP==1250
            display('LFP saved as 1250 - downsampling to 250 for save')
            swLFP = downsample(swLFP,5);
            thLFP = downsample(thLFP,5);
            sf_LFP = sf_LFP./5;

            save(scorelfppath,'thLFP','swLFP','THchannum','SWchannum','sf_LFP');
<<<<<<< HEAD
%             delete(swlfppath,thetalfppath)
=======
            delete(swlfppath,thetalfppath)
>>>>>>> origin/master
        else
            display('LFP was not saved at 1250... bug?')
            keyboard
        end
    end
end