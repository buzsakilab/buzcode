function bz_LFPFromDat(basepath)
%assumes you are in or pointed to a directory containing subdirectories for
% various recording files from a single session


%% Input and directory handling 
if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);


%% Setup
load(fullfile(basepath,[basename,'.SessionMetadata.mat']));

datname = fullfile(basepath,[basename '.dat']);
lfpname = fullfile(basepath,[basename '.lfp']);

nchannels = SessionMetadata.ExtracellEphys.Parameters.NumberOfChannels;

oldsamplerate = SessionMetadata.ExtracellEphys.Parameters.SampleRate;
newsamplerate = SessionMetadata.ExtracellEphys.Parameters.LfpSampleRate;
[N,D] = rat(newsamplerate/oldsamplerate);

%% resample one file into another
ResampleBinary(datname,lfpname,nchannels,N,D);



    
    