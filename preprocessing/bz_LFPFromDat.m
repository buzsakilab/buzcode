function bz_LFPFromDat(basepath)
%assumes you are in or pointed to a directory containing subdirectories for
% various recording files from a single session
%
%INPUT
%   Assumes presence of the following files:
%       basePath/baseName.dat
%   where basePath is a folder of the form: 
%       whateverPath/baseName/
%
%       (optional parameters files)
%       basePath/baseName.SessionMetadata.mat  (recommended)
%       basePath/baseName.xml
%
%If basePath not specified, tries the current path.
%
%
%OUTPUT
%Creates file:
%   basePath/baseName.lfp
%
%BWatson and DLevenstein 2017
%% Input and directory handling 
if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);
metadataname = fullfile(basepath,[basename,'.SessionMetadata.mat']);
xmlname = fullfile(basepath,[basename,'.xml']);

%% Setup
if exist(metadataname,'file')
    display('Getting info from SessionMetadata.mat')
    load(metadataname);
    
    nchannels = SessionMetadata.ExtracellEphys.Parameters.NumberOfChannels;
    oldsamplerate = SessionMetadata.ExtracellEphys.Parameters.SampleRate;
    newsamplerate = SessionMetadata.ExtracellEphys.Parameters.LfpSampleRate;
elseif exist(xmlname,'file')
    display(['Getting info from .xml, ',...
        'you might want to make a SessionMetadata.mat...'])
    parameters = LoadParameters(xmlname);
    
    nchannels = parameters.nChannels;
    oldsamplerate = parameters.rates.wideband;
    newsamplerate = parameters.lfpSampleRate;
else
    display(['No baseName.SessionMetadata.mat or baseName.xml file... resorting to user input,',...
        'you might want to make a SessionMetadata.mat...'])
    nchannels = input('How many channels in the file?');
    oldsamplerate = input('.dat SampleRate?');
    newsamplerate = input('Desired .lfp SampleRate?');
end

datname = fullfile(basepath,[basename '.dat']);
lfpname = fullfile(basepath,[basename '.lfp']);


[N,D] = rat(newsamplerate/oldsamplerate);

%% resample one file into another
warning(['This function uses ReampleBinary, which does NOT low-pass filter ',...
    'your LFP.  As a result it has high levels of spike contamination.  ',...
    'Consider updating this function to use ProcessResample() from NDM.'])
ResampleBinary(datname,lfpname,nchannels,N,D);

display('baseName.lfp file created! Huzzah!')

end



    
    