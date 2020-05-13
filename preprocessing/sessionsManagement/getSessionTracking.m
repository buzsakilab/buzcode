
function [tracking] = getSessionTracking(varargin)
% Tracking and concatenation over all session subfolder recordings
%
% USAGE
%
%   [tracking] = getSessionTracking(varargin)
%
% INPUTS
%   basePath       -(default: pwd) basePath for the recording file, in buzcode format:
%   roiTracking    - 2 x R, where 1C is x and 2C is y. By default it
%                   considers the whole video. With the option 'manual' allows to draw
%                   a ROI.
%   roiLED         - 2 x R. 'manual' for drawing the ROI.
%   roisPath       - provide a path with ROI mat files ('roiTRacking.mat'
%                   and 'roiLED.mat'). By default try to find it in
%                   basePath or upper folder.
%   convFact       - Spatial conversion factor (cm/px). If not provide,
%                   normalize maze size.
%   saveMat        - default true
%   forceReload    - default false
%
% OUTPUT
%       - tracking.behaviour output structure, with the fields:
%   x               - x position in cm/ normalize
%   y               - y position in cm/ normalize
%   timestamps      - in seconds, if Bazler ttl detected, sync by them
%   folder          - 
%   sync.sync       - Rx1 LED luminance.
%   sync.timestamps - 2xC with start stops of sync LED.
%
%   Manuel Valero 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'convFact',[],@isnumeric); % 0.1149
addParameter(p,'roiTracking',[],@ismatrix);
addParameter(p,'roiLED',[],@ismatrix);
addParameter(p,'roisPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'forceReload',false,@islogical)

parse(p,varargin{:});
basepath = p.Results.basepath;
convFact = p.Results.convFact;
roiTracking = p.Results.roiTracking;
roiLED = p.Results.roiLED;
roisPath = p.Results.roisPath;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) || forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end

cd(basepath); cd ..; upBasepath = pwd; cd(basepath);
if isempty(roisPath)
    if exist([basepath filesep 'roiTracking.mat'],'file') || ...
        exist([basepath filesep 'roiLED.mat'],'file')
            roisPath = basepath;
            load([roisPath filesep 'roiLED.mat'],'roiLED');
            load([roisPath filesep 'roiTracking.mat'],'roiTracking');
    elseif exist([upBasepath filesep 'roiTracking.mat'],'file') || ...
        exist([upBasepath filesep 'roiLED.mat'],'file')
            roisPath = upBasepath;
            load([roisPath filesep 'roiLED.mat'],'roiLED');
            load([roisPath filesep 'roiTracking.mat'],'roiTracking');
    end   
end

%% Find subfolder recordings
cd(basepath);
[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
C = strsplit(sessionInfo.session.name,'_');
sess = dir(strcat(C{1},'_',C{2},'*')); % get session files
count = 1;
for ii = 1:size(sess,1)
    if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Basler*avi']))
        cd([basepath filesep sess(ii).name]);
        fprintf('Computing tracking in %s folder \n',sess(ii).name);
        tempTracking{count}= LED2Tracking([],'convFact',convFact,'roiTracking',...
            roiTracking,'roiLED',roiLED,'forceReload',forceReload); % computing trajectory
        trackFolder(count) = ii; 
        count = count + 1;
    end
end
cd(basepath);

%% Concatenate and sync timestamps
ts = []; subSessions = []; maskSessions = [];
if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    for ii = 1:length(trackFolder)
        if strcmpi(MergePoints.foldernames{trackFolder(ii)},tempTracking{ii}.folder)
            sumTs = tempTracking{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1);
            subSessions = [subSessions; MergePoints.timestamps(trackFolder(ii),1:2)];
            maskSessions = [maskSessions; ones(size(sumTs))*ii];
            ts = [ts; sumTs];
        else
            error('Folders name does not match!!');
        end
    end
else
    warning('No MergePoints file found. Concatenating timestamps...');
    for ii = 1:length(trackFolder)
        sumTs = max(ts)+ tempTracking{ii}.timestamps;
        subSessions = [subSessions; [sumTs(1) sumTs(end)]];
        ts = [ts; sumTs];
    end
end

% Concatenating tracking fields...
x = []; y = []; folder = []; samplingRate = []; description = [];
for ii = 1:size(tempTracking,2) 
    x = [x; tempTracking{ii}.position.x]; 
    y = [y; tempTracking{ii}.position.y]; 
    folder{ii} = tempTracking{ii}.folder; 
    samplingRate = [samplingRate; tempTracking{ii}.samplingRate];  
    description{ii} = tempTracking{ii}.description;
end

tracking.position.x = x;
tracking.position.y = y;
tracking.folders = folder;
tracking.samplingRate = samplingRate;
tracking.timestamps = ts;
tracking.events.subSessions =  subSessions;
tracking.events.subSessionsMask = maskSessions;

[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
if saveMat
    save([basepath filesep sessionInfo.FileName '.Tracking.Behavior.mat'],'tracking');
end

end
