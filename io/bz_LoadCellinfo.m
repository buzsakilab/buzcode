function [ cellinfo,filename ] = bz_LoadCellinfo(basePath,cellinfoName,varargin)
%[ cellinfo,filename ] = bz_LoadCellinfo(basePath,cellinfoName) function for
%for loading cellinfo.mat files. cellinfo.mat files are saved as:
% datasetPath/baseName/baseName.cellinfoName.cellinfo.mat
%
%INPUT
%   basePath
%   cellinfoName    can also be 
%                   -empty, which allows use prompts the user with a list 
%                     of available cellinfo.mat files in basePath
%                   -'all' (not yet functional) load all cellinfo.mat files 
%                     for a given recording. 
% (optional inputs)
%       'UIDs'      subset of UIDs to load. (not yet functional... add this)
%       'dataset'   logical (default: false) used if basePath is a
%                   high-level dataset path. bz_LoadCellinfo then allows
%                   the user to select basepaths in the folder to load the
%                   cellinfo file from
%                   future update: 'select', 'all'
%       'baseNames' list of basepaths to load from, only used if 'dataset', true
%       'catall'    logical (default: false) if loading multiple cellinfo
%                   files from a dataset, will try to concatenate all units
%                   into a single cellinfo structure. Removes fields that
%                   are not common to all basePaths
%
%OUTPUT
%   cellinfo        loaded cellinfo structure
%   filename        filename loaded
%
%DLevenstein 2018

% Added support for NWB: Konstantinos Nasiotis 2019
%%
p = inputParser;
addParameter(p,'dataset',false)
addParameter(p,'baseNames',[])
addParameter(p,'catall',false)
parse(p,varargin{:})
dataset = p.Results.dataset;
catall = p.Results.catall;
baseNames_keep = p.Results.baseNames;


%% For loading all cellinfo files of same name from dataset

if dataset
    %Figure out which basePaths to look at
    if isempty(baseNames_keep)
        select = true;
    else
        select = false;
    end
    [basePaths,baseNames] = bz_FindBasePaths(basePath,'select',select);
    
    %Only Keep baseNames passed in 
    if ~isempty(baseNames_keep)
        keepbaseNames = ismember(baseNames,baseNames_keep);
        baseNames = baseNames(keepbaseNames);
        basePaths = basePaths(keepbaseNames);
    end
    
    %Go through each and load the cell info
    FIELDMISMATCH=false;
    for rr = 1:length(baseNames)
        thiscellinfo = bz_LoadCellinfo(basePaths{rr},cellinfoName);
        
        %Add baseName to the cellinfo file. this could be for each unit....
        thiscellinfo.baseName = repmat(baseNames(rr),size(thiscellinfo.UID));
        
        %Check if the new .mat has any additional fields
        if exist('cellinfo','var')    
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
            [cellinfo,thiscellinfo] = bz_Matchfields(cellinfo,thiscellinfo,'remove');
        end
        
        cellinfo(rr) = thiscellinfo;
    end
    
    if FIELDMISMATCH
        warning('One or more of your .mats has missing fields')
    end
    
    if catall
        cellinfo = bz_CollapseStruct(cellinfo,'match','justcat',true);
    end
    
    return %send out the compiled cellinfo structure
end


%% For loading cellinfo from a single basePath
if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);

if ~exist('cellinfoName','var') || isempty(cellinfoName)
    allCellinfoFiles = dir(fullfile(basePath,[baseName,'.','*','.cellinfo.mat']));
    
    if ~isempty(allCellinfoFiles)
        [s,v] = listdlg('PromptString','Which cellinfo.mat would you like to load?',...
                     'ListString',{allCellinfoFiles.name},'SelectionMode','single');
        if isempty(s)
            cellinfo = []; filename = [];
            return
        end
        filename = fullfile(basePath,allCellinfoFiles(s).name);
    else
        filename = [];
    end

else
    filename = fullfile(basePath,[baseName,'.',cellinfoName,'.cellinfo.mat']);
end


if exist(filename,'file')
    cellinfostruct = load(filename);
else % Check for NWB file
    
    cellinfostruct = loadCellInfoNWB(basePath);
    
    if isempty(cellinfostruct)
        cellinfo = [];
        return
    end
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




%% Function that loads NWB fields to a Behavior structure
function cellinfostruct = loadCellInfoNWB(basePath)
    
    % Check for NWB file since no .behavior.mat files were found
    d = dir(fullfile(basePath, '*.nwb'));
    if length(d) > 1 % If more than one .nwb files exist in the directory, select which one to load from
        warning('there is more than one .nwb file in this directory');
        [iNWBFile, ~] = listdlg('PromptString','Which NWB file would you like to load?',...
                                 'ListString',{d.name},'SelectionMode','single');
        d = d(iNWBFile);
    elseif length(d) == 0
        d = dir('*nwb');
        if isempty(d)
            disp('No .cellinfo.mat or NWB files present in this directory')
            cellinfostruct = [];
            return
        end
    end
    
    % If an NWB file exists, load the info according the Buzcode standard
    nwb_file = fullfile(d.folder, d.name); 

    nwb2 = nwbRead(nwb_file);
    
    if ~isempty(nwb2.units)
    
        nNeurons = length(nwb2.units.id.data.load);

        cellinfostruct.spikes.UID           = double(nwb2.units.id.data.load');
        cellinfostruct.spikes.shankID       = nwb2.units.vectordata.get('shankID').data.load';
        cellinfostruct.spikes.cluID         = nwb2.units.vectordata.get('cluID').data.load';
        cellinfostruct.spikes.maxWaveformCh = nwb2.units.vectordata.get('maxWaveformCh').data.load';
        cellinfostruct.spikes.region        = nwb2.units.vectordata.get('region').data.load';
        cellinfostruct.spikes.sessionName   = nwb2.identifier;

        for iNeuron = 1:nNeurons
            cellinfostruct.spikes.rawWaveform{1,iNeuron} = nwb2.units.waveform_mean.data.load([iNeuron, 1], [iNeuron,Inf]);

            if iNeuron == 1 
                cellinfostruct.spikes.times{1,iNeuron} = nwb2.units.spike_times.data.load(1, double(nwb2.units.spike_times_index.data.load(iNeuron)));
            else
                cellinfostruct.spikes.times{1,iNeuron} = nwb2.units.spike_times.data.load(double(nwb2.units.spike_times_index.data.load(iNeuron-1))+1, double(nwb2.units.spike_times_index.data.load(iNeuron)));
            end
        end

        %% Optional
        % reorder for easy comparison with the tutorial example
        C = {'UID','times','shankID','cluID','rawWaveform','maxWaveformCh','region','sessionName'};
        cellinfostruct.spikes = orderfields(cellinfostruct.spikes,C);
    
    else
        disp('No units information in the selected NWB file');
        cellinfostruct = [];
        return
    end

end














end


