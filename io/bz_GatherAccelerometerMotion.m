function motion = bz_GetAccelerometerMotion(basepath,AccelerChans,varargin)
%function motion = bz_GetAccelerometerMotion(basepath,AccelerChans,varargin)
%
% Grabs data from Intan accelerometers and outputs normalized movement per
%time bin.  Time bin is specified by OutputSampFreq in seconds, with a
%default of 2 meaning 2Hz or 500ms bins.  
%
%
% INPUT
% Required inputs:
%   basepath        folder containing .xml and .lfp files.
%                   basepath and files should be of the form:
%                   'whateverfolder/recordingName/recordingName'
%   'AccelerChans'  Channels devoted to motion/accelerometer.  Needed only
%                   if the option 'MotionSource" has the value "Accelerometer' 
%
% Optional inputs: 
%   'OutputSampFreq' Frequency in Hz of bins for output.  Each bin will
%                   have averaged movement per bin.  Default is
%                   OutputSampFreq=2 meaning 2 bins per second
%   'savedir'       Default: basepath
%   'savebool'      Default: true. Save anything.
%
% OUTPUT
%   motion          A struct array containing: 
%           - .data - vector of motion values per bin
%           - .timestamps - timestamp of middle-point of bin, in seconds
%               from start of recording
%           - .data_var: a version of teh data based on moving
%               variance, sometimes more useful than raw movement
%           - .AccelerChans - channels used and averaged for accelerometer
%               data
%           - .detectorName - 'IntanAccelerometer' default for now


%% Handle case of no input
if ~exist('basepath','var')
    basepath = uigetdir(cd,...
        'Which recording(s) would you like to state score?');
    if isequal(basepath,0);return;end  
end
if ~exist('AccelerChans','var')
    error('Must input accelerometer channel numbers for them to be loaded')
end

basename = bz_BasenameFromBasepath(basepath);

%% Parse inputs
p = inputParser;
addParameter(p,'OutputSampFreq',2,@isnumeric)
addParameter(p,'savebool',true,@islogical)
addParameter(p,'savedir',basepath)

parse(p,varargin{:})
OutputSampFreq = p.Results.OutputSampFreq;
savedir = p.Results.savedir;
savebool = p.Results.savebool;

%% Finish setting up
savefolder = fullfile(savedir,basename);
if ~exist(savefolder,'dir')
    mkdir(savefolder)
end
sessionInfo = bz_getSessionInfo(basepath,'noPrompts',true);
lfppath = fullfile(basepath,[sessionInfo.session.name,'.lfp']);
savepath = fullfile(savefolder,[basename '.AccelerometerMotion.mat']);


%% Load data and calc normalized total motion
motiondata = bz_LoadBinary(lfppath, 'channels', AccelerChans+1, 'nChannels', sessionInfo.nChannels)';

motiondata = single(motiondata);        
motiondata = abs(zscore(motiondata')');%total relative movement per channel
motiondata = sum(motiondata, 1);%summate across channels

%% "Binning"
Fs = sessionInfo.lfpSampleRate;
binScootS = 1 ./ OutputSampFreq;
binScootSamps = round(Fs*binScootS); % must be integer, or error             
%         bin_numsamps = round(binScootS*Fs);
%         bin_inds = -bin_numsamps:bin_numsamps;%+- that number of ms in samples
%         timestamps = (1+bin_inds(end)):binScootSamps:(size(motiondata,2)-bin_inds(end));
%         numbins = length(timestamps);

%mean in each bin
newlen = size(motiondata,2) - mod(size(motiondata,2),binScootSamps);
motiondata = motiondata(:,1:newlen);   
motiondata = reshape(motiondata,[binScootSamps size(motiondata,2)/binScootSamps]);
motiondata = mean(motiondata,1);

timestamps = 1:newlen;
timestamps = reshape(timestamps,[binScootSamps size(timestamps,2)/binScootSamps]);
timestamps = mean(timestamps,1)/1250;

%take the variance over each chunk of 3 seconds
motiondatavar = movvar(motiondata,3);%show the variance across 3 points, not the raw value. 
    %movvar function slides one sample at a time, not a full bin width

%% Function output
motion.data = motiondata;
motion.timestamps = timestamps;
motion.data_var = motiondatavar;
motion.channels = AccelerChans;
motion.detectorName = 'IntanAccelerometer';
motion.OutputSampFreq = OutputSampFreq;

%% Save to disk 
if savebool
    save(savepath,'motion')
end
