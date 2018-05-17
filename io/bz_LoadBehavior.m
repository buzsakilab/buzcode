function [ behavior,filename ] = bz_LoadBehavior( basePath,behaviorName )
%[ behavior,filename ] = bz_LoadBehavior( basePath,behaviorName ) function for
%loading behavior.mat files.
%
%
%DLevenstein 2017
%%

if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);

if ~exist('behaviorName','var')
    allBehaviorFiles = dir(fullfile(basePath,[baseName,'.','*','.behavior.mat']));
    [s,v] = listdlg('PromptString','Which behavior.mat would you like to load?',...
                 'ListString',{allBehaviorFiles.name},'SelectionMode','single');
    if isempty(s)
        behavior = []; filename = [];
        return
    end
    filename = fullfile(basePath,allBehaviorFiles(s).name);
else
    filename = fullfile(basePath,[baseName,'.',behaviorName,'.behavior.mat']);
end


if exist(filename,'file')
    behaviorstruct = load(filename);
else
    warning([filename,' does not exist...'])
    behavior = [];
    return
end

varsInFile = fieldnames(behaviorstruct);

if numel(varsInFile)==1
    behavior = behaviorstruct.(varsInFile{1});
else
    warning('Your .behavior.mat has multiple variables/structures in it... wtf.')
    behavior = behaviorstruct;
end

%Check that the behavior structure meets buzcode standards
[isBehavior] = bz_isBehavior(behavior);
switch isBehavior
    case false
        warning('Your behavior structure does not meet buzcode standards. Sad.')
end


end

