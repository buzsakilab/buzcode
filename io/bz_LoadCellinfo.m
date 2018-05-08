function [ cellinfo,filename ] = bz_LoadCellinfo(basePath,cellinfoName)
%[ events ] = bz_LoadCellinfo(basePath,eventsName) function for
%loading cellinfo.mat files. cellinfo.mat files are saved as...
% datasetPath/baseName/baseName.cellinfoName.cellinfo.mat
%
%cellinfoName can be the name of a cellinfo.mat file, If empty, prompts the user
%with a list of available cellinfo.mat files in basePath.
%Future update: 'all' (nonfunctional) to load all cellinfo.mat files for a given recording. 
%
%DLevenstein 2018
%%
if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);

if ~exist('cellinfoName','var')
    allCellinfoFiles = dir(fullfile(basePath,[baseName,'.','*','.cellinfo.mat']));
    [s,v] = listdlg('PromptString','Which cellinfo.mat would you like to load?',...
                 'ListString',{allCellinfoFiles.name},'SelectionMode','single');
    if isempty(s)
        cellinfo = []; filename = [];
        return
    end
    filename = fullfile(basePath,allCellinfoFiles(s).name);
else
    filename = fullfile(basePath,[baseName,'.',cellinfoName,'.cellinfo.mat']);
end


if exist(filename,'file')
    cellinfostruct = load(filename);
else
    warning([filename,' does not exist...'])
    cellinfo = [];
    return
end

varsInFile = fieldnames(cellinfostruct);

if numel(varsInFile)==1
    cellinfo = cellinfostruct.(varsInFile{1});
else
    warning('Your .cellinfo.mat has multiple variables/structures in it... wtf.')
    cellinfo = cellinfostruct;
end

%Check that the events structure meets buzcode standards
[isCellinfo] = bz_isCellInfo(cellinfo);
switch isCellinfo
    case false
        warning('Your cellinfo structure does not meet buzcode standards. Sad.')
end

end


