
function [sessionArmChoice] = getSessionArmChoice(varargin)
% Compute and get session arm choice over all session subfolders
%
% USAGE
%
%   [tracking] = getSessionTracking(varargin)
%
% INPUTS
% basePath                      (default: pwd) basePath for the recording file, 
%                                    in buzcode format:
% task                          'alternation' and 'cudeSide'
% forceReload                   Force detection (boolean, default false)
% verbose                       Default false
% saveMat                       Default true
% forceRun                      Try to run armChoice even if there is no
%                                   tracking files on the directories (default false).
%
% OUTPUT
%       - armChoice.behaviour.(subSessionFolder) output structure, with the fields:
% armChoice.timestamps          Choice timestamps, in seconds
% armChoice.arm                 Choosed arm, 0 is left, 1 is right
% armChoice.delay.ints          Delay intervals, in seconds
% armChoice.delay.dur           Delay duration, in seconds
% armChoice.delay.timestamps    Delay timestamps, in seconds
% armChoice.choice              Performance vector, 1 is right choice, 0 is
%                                   wrong. First choice is Nan.
% armChoice.performance         Alternation probability (#alternation/#trials)
% armChoice.forzed              1 if forzed alternation, 0 if spontaneous
%                                   alternation
% armChoice.task                'alternation' and 'cudeSide'  
%
%   Manuel Valero 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'task',[],@ischar);
addParameter(p,'force',false,@islogical);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'forceRun',false,@islogical)
parse(p,varargin{:});
task = p.Results.task;
forceReload = p.Results.force;
basepath = p.Results.basepath;
verbose = p.Results.verbose;
saveMat = p.Results.saveMat;
forceRun = p.Results.forceRun;


%% Deal with inputs
if ~isempty(dir([basepath filesep '*.SessionArmChoice.Events.mat'])) || forceReload
    disp('Session arm choice already detected! Loading file.');
    file = dir([basepath filesep '*.SessionArmChoice.Events.mat']);
    load(file.name);
    return
end

%% Find subfolder recordings
cd(basepath);
try [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
    C = strsplit(sessionInfo.session.name,'_');
    sess = dir(strcat(C{1},'_',C{2},'*')); % get session files
catch
    warning('No sessionInfo found!');
    sess = dir(pwd);
    sess(1:2) = [];
end

count = 1;
for ii = 1:size(sess,1)
    if sess(ii).isdir
        if ~isempty(dir([basepath filesep sess(ii).name filesep '*Basler*avi'])) || forceRun 
            cd([basepath filesep sess(ii).name]);
            fprintf('Computing arm Choice in %s folder \n',sess(ii).name);
            sessionArmChoice.(sess(ii).name)= getArmChoice('verbose',verbose,'task',task);
            trackFolder(count) = ii; 
            count = count + 1;
        end
    end
end
cd(basepath);

efields = fieldnames(sessionArmChoice);
try tracking = getSessionTracking;
    if size(tracking.events.subSessions,1) == size(efields,1)
    disp('Correctiong timestamps for session recording...');
    for ii = 1:size(efields,1)
        preRec = tracking.events.subSessions(ii,1);
        sessionArmChoice.(efields{ii}).timestamps = ...
            sessionArmChoice.(efields{ii}).timestamps + preRec;
        sessionArmChoice.(efields{ii}).delay.timestamps = ...
            sessionArmChoice.(efields{ii}).delay.timestamps + preRec;
    end
else
    warning('Number of behavioral recordings do not match!')
end
catch 
    warning('No available tracking!');
end

try [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
catch
    sessionInfo.FileName = split(pwd,filesep); sessionInfo.FileName = sessionInfo.FileName{end};
end
if saveMat
    save([basepath filesep sessionInfo.FileName '.SessionArmChoice.Events.mat'],'sessionArmChoice');
end

end

