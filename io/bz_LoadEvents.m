function [ events ] = bz_LoadEvents(basePath,eventsName)
%[ events ] = bz_LoadEvents(basePath,eventsName) function for
%loading events.mat files. events.mat files are saved as...
% datasetPath/baseName/baseName.eventsName.events.mat
%
%eventsName can be the name of a events.mat file, or can be 'all' to load
%all events.mat files for a given recording. If empty, prompts the user
%with a list of available events.mat files in basePath
%
%DLevenstein 2017
%%
baseName = bz_BasenameFromBasepath(basePath);

%if strcmp('eventsName','all')
%     allEventsFiles = dir(fullfile(basePath,[baseName,'.','*','.events.mat']));
%         [s,v] = listdlg('PromptString','Which events.mats would you like to load?',...
%                         'ListString',allEventsFiles);
%         eventsfilestoload = allEventsFiles(s);
        
%end


eventsfile = fullfile(basePath,[baseName,'.',eventsName,'.events.mat']);


%evalin('caller',['load(''',eventsfile,''')']);

if exist(eventsfile,'file')
    eventstruct = load(eventsfile);
else
    warning([eventsfile,' does not exist...'])
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

end

