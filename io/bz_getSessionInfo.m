function [sessionInfo] = bz_getSessionInfo(basePath,varargin)

%% inputs and defaults
p = inputParser;
addParameter(p,'noPrompts',false,@islogical);
parse(p,varargin{:})
noPrompts = p.Results.noPrompts;

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
if ~isfield(sessionInfo,'region')
    regionadd = questdlg(['Your sessionInfo is missing regions, ',...
        'would you like to add them?']);
    if strcmp(regionadd,'Yes')
        [ sessionInfo ] = bz_sessionInfoGUI( sessionInfo );
        SIexist = false; 
    end
end
    
    
%% Save sessionInfo file   
%Here: prompt user to save basePath/baseName.sessionInfo.mat if loaded from
%xml
if ~noPrompts && ~SIexist %Inform the user that they should save a file for later
    savebutton = questdlg(['Would you like to save your sessionInfo in ',...
        filename]);
    if strcmp(savebutton,'Yes'); save(filename,'sessionInfo'); end
end

end