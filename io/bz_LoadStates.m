function [ states,filename ] = bz_LoadStates(basePath,statesName)
%[ states ] = bz_LoadStates(basePath,statesName) function for
%loading states.mat files. states.mat files are saved as...
% datasetPath/baseName/baseName.statesName.states.mat
%
%statesName can be the name of a states.mat file, or can be 'all' to load
%all states.mat files for a given recording. If empty, prompts the user
%with a list of available states.mat files in basePath
%
%DLevenstein 2017
%%
baseName = bz_BasenameFromBasepath(basePath);

%Here is where the code will go to allow user to select states or load all
%states
%if strcmp('statesName','all')
%     allStatesFiles = dir(fullfile(basePath,[baseName,'.','*','.states.mat']));
%         [s,v] = listdlg('PromptString','Which states.mats would you like to load?',...
%                         'ListString',allStatesFiles);
%         statesfilestoload = allStatesFiles(s);
        
%end


filename = fullfile(basePath,[baseName,'.',statesName,'.states.mat']);

if exist(filename,'file')
    statestruct = load(filename);
else
    warning([filename,' does not exist...'])
    states = [];
    return
end

varsInFile = fieldnames(statestruct);

if numel(varsInFile)==1
    states = statestruct.(varsInFile{1});
else
    warning('Your .states.mat has multiple variables/structures in it... wtf.')
    states = statestruct;
end

%Check here that the structure complies with buzcode format?
end

