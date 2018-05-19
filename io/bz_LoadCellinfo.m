function [ cellinfo,filename ] = bz_LoadCellinfo(basePath,cellinfoName,varargin)
%[ events ] = bz_LoadCellinfo(basePath,eventsName) function for
%loading cellinfo.mat files. cellinfo.mat files are saved as...
% datasetPath/baseName/baseName.cellinfoName.cellinfo.mat
%
%cellinfoName can be the name of a cellinfo.mat file, If empty, prompts the user
%with a list of available cellinfo.mat files in basePath.
%Future update: 'all' (nonfunctional) to load all cellinfo.mat files for a given recording. 
%
% (optional inputs)
%       'dataset'   logical (default: false) used if basePath is a
%                   high-level dataset path. bz_LoadCellinfo then allows
%                   the user to select basepaths in the folder to load the
%                   cellinfo file from
%                   future update: 'select', 'all'
%
%DLevenstein 2018
%%
p = inputParser;
addParameter(p,'dataset',false,@islogical)
parse(p,varargin{:})
dataset = p.Results.dataset;

%% For loading all cellinfo files of same name from dataset

if dataset
    %Figure out which basePaths to look at
    [basePaths,baseNames] = bz_FindBasePaths(basePath,...
        'select',true);
    
    %Go through each and load the cell info
    FIELDMISMATCH=false;
    numrecs = length(baseNames);
    for rr = 1:numrecs
        thiscellinfo = bz_LoadCellinfo(basePaths{rr},cellinfoName);
        
        %Add baseName to the cellinfo file
        thiscellinfo.baseName = baseNames{rr};
        
        %Check if the new .mat has any additional fields
        if exist('cellinfo','var')    
            matfields = fieldnames(thiscellinfo);
            resultsfields = fieldnames(cellinfo);
            newfields = setdiff(matfields,resultsfields);
            if ~isempty(newfields)
                for ff = 1:length(newfields)
                    cellinfo(1).(newfields{ff}) = []; 
                end
                FIELDMISMATCH=true;
            end

            cellinfo = orderfields(cellinfo,thiscellinfo);   
        end
        
        cellinfo(rr) = thiscellinfo;
    end
    
    if FIELDMISMATCH
        disp('One or more of your .mats has missing fields')
    end
    
    %Concatenate structures here. find matching dimensions and cat along
    %the non-matching?
    return %send out the compiled cellinfo structure
end


%% For loading cellinfo from a single basePath
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


