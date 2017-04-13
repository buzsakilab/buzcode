function bz_LFPFromDat(basepath,basename)
%assumes you are in or pointed to a directory containing subdirectories for
% various recording files from a single session


%% Input and directory handling 
if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end

if ~exist('basename','var')
    [~,basename] = fileparts(basepath);
elseif isempty(basename)
    [~,basename] = fileparts(basepath);
end

%% Setup
load(fullfile(basepath,[basename,'_SessionMetadata.mat']));

datname = fullfile(basepath,[basename '.dat']);
lfpname = fullfile(basepath,[basename '.lfp']);

nchannels = SessionMetadata.ExtracellEphys.Parameters.NumberOfChannels;

oldsamplerate = SessionMetadata.ExtracellEphys.Parameters.SampleRate;
newsamplerate = SessionMetadata.ExtracellEphys.Parameters.LfpSampleRate;
[N,D] = rat(newsamplerate/oldsamplerate);

%% resample one file into another
ResampleBinary(datname,lfpname,nchannels,N,D);



    
    