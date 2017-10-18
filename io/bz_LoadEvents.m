function [ events,filename ] = bz_LoadEvents(basePath,eventsName)
%[ events ] = bz_LoadEvents(basePath,eventsName) function for
%loading events.mat files. events.mat files are saved as...
% datasetPath/baseName/baseName.eventsName.events.mat
%
%eventsName can be the name of a events.mat file, If empty, prompts the user
%with a list of available events.mat files in basePath.
%Future update: 'all' (nonfunctional) to load all events.mat files for a given recording. 
%
%DLevenstein 2017
%%
if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);

if ~exist('eventsName','var')
    allEventsFiles = dir(fullfile(basePath,[baseName,'.','*','.events.mat']));
    [s,v] = listdlg('PromptString','Which events.mat would you like to load?',...
                 'ListString',{allEventsFiles.name},'SelectionMode','single');
    if isempty(s)
        events = []; filename = [];
        return
    end
    filename = fullfile(basePath,allEventsFiles(s).name);
else
    filename = fullfile(basePath,[baseName,'.',eventsName,'.events.mat']);
end


if exist(filename,'file')
    eventstruct = load(filename);
else
    warning([filename,' does not exist...'])
    events = [];
    return
end

varsInFile = fieldnames(eventstruct);

if numel(varsInFile)==1
    events = eventstruct.(varsInFile{1});
else
    warning('Your .events.mat has multiple variables/structures in it... wtf.')
    events = eventstruct;
end

%Check that the events structure meets buzcode standards
[isEvents] = bz_isEvents(events);
switch isEvents
    case false
        warning('Your events structure does not meet buzcode standards. Sad.')
end

end

