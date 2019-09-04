function [ analysisResults,filename ] = bz_LoadAnalysisResults(basePath,analysisName,varargin)
%[ cellinfo,filename ] = bz_LoadCellinfo(basePath,cellinfoName) function for
%for loading cellinfo.mat files. cellinfo.mat files are saved as:
% datasetPath/baseName/baseName.cellinfoName.cellinfo.mat
%
%INPUT
%   basePath
%   analysisName    can also be 
%                   -empty, which allows use prompts the user with a list 
%                     of available AnalysisResults.XXXXXXX.mat files in basePath
%                   -'all' (not yet functional) load all cellinfo.mat files 
%                     for a given recording. 
% (optional inputs)
%       'UIDs'      subset of UIDs to load. (not yet functional... add this)
%       'dataset'   logical (default: false) used if basePath is a
%                   high-level dataset path. bz_LoadAnalysisResults then allows
%                   the user to select basepaths in the folder to load the
%                   cellinfo file from
%                   future update: 'select', 'all'
%       'baseNames' list of basepaths to load from, only used if 'dataset', true
%                   use 'all' to load from all basePaths without prompt
%       'catall'    logical (default: false) if loading multiple cellinfo
%                   files from a dataset, will try to concatenate all units
%                   into a single cellinfo structure. Removes fields that
%                   are not common to all basePaths
%
%OUTPUT
%   analysisResults         loaded analysisResults structure
%   filename                filename loaded
%
%DLevenstein 2018
%%
p = inputParser;
addParameter(p,'dataset',false)
addParameter(p,'baseNames',[])
addParameter(p,'catall',false)
parse(p,varargin{:})
dataset = p.Results.dataset;
catall = p.Results.catall;
baseNames_keep = p.Results.baseNames;


%% For loading all AnalysisResults files of same name from dataset

if dataset
    %Figure out which basePaths to look at
    if isempty(baseNames_keep)
        select = true;
    else
        select = false;
    end
    [basePaths,baseNames] = bz_FindBasePaths(basePath,'select',select);
    
    %Only Keep baseNames passed in 
    if ~isempty(baseNames_keep) && ~strcmp(baseNames_keep,'all')
        keepbaseNames = ismember(baseNames,baseNames_keep);
        baseNames = baseNames(keepbaseNames);
        basePaths = basePaths(keepbaseNames);
    end
    
    %Go through each and load the cell info
    FIELDMISMATCH=false;
    for rr = 1:length(baseNames)
        thisAResults = bz_LoadAnalysisResults(basePaths{rr},analysisName);
        
        %Add baseName to the cellinfo file. this could be for each unit....
        thisAResults.baseName = baseNames(rr);
        
        %Check if the new .mat has any additional fields
        if exist('analysisResults','var')    
%             matfields = fieldnames(thiscellinfo);
%             resultsfields = fieldnames(cellinfo);
%             newfields = setdiff(matfields,resultsfields);
%             if ~isempty(newfields)
%                 for ff = 1:length(newfields)
%                     cellinfo(1).(newfields{ff}) = []; 
%                 end
%                 FIELDMISMATCH=true;
%             end
% 
%             cellinfo = orderfields(cellinfo,thiscellinfo);   
            [analysisResults,thisAResults] = bz_Matchfields(analysisResults,thisAResults,'remove');
        end
        
        analysisResults(rr) = thisAResults;
    end
    
    if FIELDMISMATCH
        warning('One or more of your .mats has missing fields')
    end
    
    if catall
        analysisResults = bz_CollapseStruct(analysisResults,'match','justcat',true);
    end
    
    filename = baseNames;
    return %send out the compiled cellinfo structure
end


%% For loading cellinfo from a single basePath
if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);

if ~exist('analysisName','var') || isempty(analysisName)
    allAResultsFiles = dir(fullfile(basePath,[baseName,'.AnalysisResults.','*','.mat']));
    [s,v] = listdlg('PromptString','Which AnalysisResults would you like to load?',...
                 'ListString',{allAResultsFiles.name},'SelectionMode','single');
    if isempty(s)
        analysisResults = []; filename = [];
        return
    end
    filename = fullfile(basePath,allAResultsFiles(s).name);
else
    filename = fullfile(basePath,[baseName,'.AnalysisResults.',analysisName,'.mat']);
end


if exist(filename,'file')
    analysisResults = load(filename);
else
    warning([filename,' does not exist...'])
    analysisResults = [];
    return
end


end


