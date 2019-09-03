function [sessionInfo] = bz_getSessionInfo(basePath,varargin)
%[sessionInfo] = bz_getSessionInfo(basePath) loads the sessionInfo metadata
%for the recording in basePath. basePath should be in the format:
%       /whateverPath/baseName/
%           a file  basePath/baseName.sessionInfo.mat
%           or      basePath/baseName.xml
%           should exist.
%If no baseName.sessionInfo.mat exists, loads from the xml.
%
%INPUT
%   basePath            directory: '/whatevetPath/baseName/'
%   (options)
%       'saveMat'       (default: prompt)
%       'noPrompts'     (default: false) prevents prompts about
%                       saving/adding metadata
%       'editGUI'       (default: false) opens a GUI to edit select
%                       sessionInfo fields (beta, please add/improve!)
%
%OUTPUT
%   sessionInfo         metadata structure
%
%2017 DLevenstein and DTingley
%% inputs and defaults
p = inputParser;
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'editGUI',false,@islogical);
addParameter(p,'saveMat',false,@islogical);
parse(p,varargin{:})
noPrompts = p.Results.noPrompts;
editGUI = p.Results.editGUI;
saveMat = p.Results.saveMat;

if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);
filename = fullfile(basePath,[baseName,'.sessionInfo.mat']);

%% Load the stuff
%d = dir('*sessionInfo*'); %all files with sessioninfo in the name
if exist(filename,'file')
    sessionInfostruct = load(filename);
    %Checks that there is a single structure in the sessionInfo file
    varsInFile = fieldnames(sessionInfostruct); 
    if numel(varsInFile)==1
        sessionInfo = sessionInfostruct.(varsInFile{1});
    else
        warning('Your .sessionInfo.mat has multiple variables/structures in it... wtf.')
        sessionInfo = sessionInfostruct;
    end
    SIexist = true;  %Marks that session info exists as expected
else
   warning(['could not find file ',baseName,'.sessionInfo.mat ',...
       'running LoadParameters instead..']) 
   sessionInfo = LoadParameters(basePath);
   SIexist = false; 
end

%% Check sessionInfo using bz_isSessionInfo and update if necessary
bz_isSessionInfo(sessionInfo);

%Here: prompt user to add any missing sessionInfo fields and save
if editGUI
    [ sessionInfo ] = bz_sessionInfoGUI(sessionInfo);
    SIexist = false;
% elseif ~isfield(sessionInfo,'region') && ~noPrompts
%     regionadd = questdlg(['Your sessionInfo is missing regions, ',...
%         'would you like to add them?'],'Add Regions?','Yes');
%     switch regionadd
%         case 'Yes'
%             [sessionInfo] = bz_sessionInfoGUI(sessionInfo,'Regions');
%             SIexist = false; 
%         case 'Cancel'
%             return
%     end
end

%Should check that sessionInfo.session.name and sesioninfo.session.path
%match basePath....  if not, prompt the user to save with the correct
%information.
    
%% Save sessionInfo file   
%Prompt user to save basePath/baseName.sessionInfo.mat 
%if loaded from xml or changed
if ~noPrompts && ~SIexist %Inform the user that they should save a file for later
    savebutton = questdlg(['Would you like to save your sessionInfo in ',...
        filename '?'],'Save sessionInfo?','Yes');
    switch savebutton
        case 'Yes'
            saveMat = true;
        case 'No'
            saveMat = false;
        case 'Cancel'
            return
    end
end

if saveMat
    disp(['saving ',baseName,'.sessionInfo.mat'])
    save(filename,'sessionInfo'); 
end
end
