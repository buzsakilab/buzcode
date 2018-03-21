function bz_LFPFromAbf(basePath,chanNums)
%assumes you are in or pointed to a directory containing subdirectories for
% various recording files from a single session
%
%INPUT
%   Assumes presence of the following files:
%       basePath/baseName.abf
%   where basePath is a folder of the form: 
%       whateverPath/baseName/
%
%       (optional parameters files)
%       basePath/baseName.SessionInfo.mat  (recommended)
%       basePath/baseName.xml
%
%   (optional)
%       chanNums    (1-indexed) channel numbers of the lfp. (will be
%                   converted to 0-indexed for ChanID)
%
%If basePath not specified, tries the current working directory.
%
%
%OUTPUT
%Creates file:
%   basePath/baseName.lfp
%
%DLevenstein 2018
%% Input and directory handling 
if ~exist('basePath','var')
    basePath = cd;
elseif isempty(basePath)
    basePath = cd;
end
baseName = bz_BasenameFromBasepath(basePath);

abffilename = fullfile(basePath,[baseName,'.abf']);

if ~exist(abffilename,'file')
    basepathfiles = dir([basePath,filesep,'*.abf']);
    if length(basepathfiles)==1
        button = questdlg(['There is no file ',abffilename,...
            ', but I did find ',basepathfiles.name,...
            '. Would you like to open that?']);
        switch button
            case 'Yes'
                abffilename = fullfile(basePath,basepathfiles.name);
            case {'No','Cancel'}
                return
        end
    elseif length(basepathfiles)>1
        display(['Multiple .abf files in basePath... try following buzcode standards, eh?',...
            '  Name your file basePath/baseName.abf'])
        return
    else
        display('No .abf file in basePath... try again later')
        return
    end
end

%% Load the abf file
[raw_data,si,file_info]=abfload(abffilename);


%% Prompt the user, which channels are LFP
numabfchans = length(raw_data(1,:));
if ~exist('chanNums','var')
    chanNums = listdlg('ListString',file_info.recChNames,...
        'PromptString',['Which channels are LFP? ']);
    %Replace this with prompt file_info.recChNames
end

%% Load the sessionInfo file
try
    sessionInfo = bz_getSessionInfo(basePath);
catch
    display('no sessionInfo file found, making one!')
    sessionInfo.FileName = baseName;
end

%Fill in sessionInfo from abf file and save
sessionInfo.nChannels = length(chanNums);
sessionInfo.chanNames = file_info.recChNames(chanNums);
sessionInfo.channels = 0:sessionInfo.nChannels - 1;
sessionInfo.Date = file_info.uFileStartDate;
sessionInfo.abfSampleRate = 1./(si*10^-6);
sessionInfo.lfpSampleRate = 1250;
sessionInfo.Amplification = 1000;
sessionInfo.spikeGroups.groups{1} = sessionInfo.channels;
sessionInfo.spikeGroups.nGroups = 1;
sessionInfo.file_info = file_info;

%Calculate the resampling ratios
[up,down] = rat(sessionInfo.lfpSampleRate/sessionInfo.abfSampleRate);

%Sweeps
numsweeps = length(file_info.sweepStartInPts);
if numsweeps>1
    raw_data = reshape(permute(raw_data,[1 3 2]),[],4);
    sessionInfo.sweeps_dt = [ceil(file_info.sweepStartInPts.*up./down) ...
        ceil((file_info.sweepStartInPts+file_info.sweepLengthInPts).*up./down)-1];
    sessionInfo.sweeps_s = (sessionInfo.sweeps_dt-1)./sessionInfo.lfpSampleRate;
end

sessionInfoname = fullfile(basePath,[baseName,'.sessionInfo.mat']);
save(sessionInfoname,'sessionInfo'); 

%% Resample and save the .lfp file
lfpname = fullfile(basePath,[baseName '.lfp']);
outputFile = fopen(lfpname,'w');

resampled = resample(raw_data(:,chanNums),up,down);
resampled = resampled.*sessionInfo.Amplification;

count = fwrite(outputFile,resampled','int16');

end

%%
