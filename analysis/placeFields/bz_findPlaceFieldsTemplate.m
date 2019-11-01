function [placeFieldTemplate] = bz_findPlaceFieldsTemplate(varargin)
%   [placeFieldTemplate] = bz_findPlaceFields1D(firingMaps)
%   Creates a template with all cells acording to the arragement of their
%   place fields in a 1D maze. 
%
%   INPUTS
%
%   placeFieldStats cellinfo structure from bz_findPlaceFields1D with the 
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%
%    'placeFieldStats'  cellinfo structure from bz_findPlaceFields1D with
%                       the following fields
%		.UID            unit ids
%		.sessionName    name of session
%		.params         parameters
%       .mapStats       Statistics of the Firing Map
%                       .x      Bin position of maximum firing rate (center
%                               of place field)
%     'firingMaps'  cellinfo struct with the following fields
%        .rateMaps              gaussian filtered rates
%        .rateMaps_unsmooth     raw rate data
%        .rateMaps_box          box filtered rates
%        .countMaps             raw spike count data
%        .occuMaps              position occupancy data
%                   ouput structure from bz_firingMapAvg. If not provided,
%                   it loads it from 'basepath' or current folder
%     'UIDs'        A Mx1 boolean matrix with 1s for units to be considered
%                   and 0s for units to be discarded.(M: number of units)
%     'basepath'    full path where session is located (default pwd)
%                   e.g. /mnt/Data/buddy140_060813_reo/buddy140_060813_reo
%     'saveMat'   	Saves file, logical (default: true) 
%    =========================================================================
%
%   OUTPUTS
%
%   placeFieldTemplate structure with the following fields:
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%
%    .Peak             1 x (# conditions) cell array. 
%                      Within each cell there is a (# units) x 3 matrix.
%                      First column has the bin corresponding to the firing 
%                      map peak for each unit (NaN for units without place
%                      field). Second column has the unit ID. Third column
%                      has the position of the unit in the UID vector. The
%                      third dimension contains this similar matrix for the
%                      other conditions.
%    .CenterofMass     1 x (# conditions) cell array. 
%                      Within each cell there is a (# units) x 3 matrix.
%                      First column has the center of mass of the firing
%                      rate map for each unit (NaN for units without place
%                      field). Second column has the unit ID. Third column
%                      has the position of the unit in the UID vector. The
%                      third dimension contains this similar matrix for the
%                      other conditions.
%    =========================================================================
%    
% Andrea Navas-Olive, Antonio FR, 10/2019

% Parse inputs 
p=inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'firingMapsAvg',{},@isstruct);
addParameter(p,'placeFieldStats',{},@isstruct);
addParameter(p,'saveMat', true, @islogical);
addParameter(p,'UIDs',[],@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
firingMaps = p.Results.firingMapsAvg;
placeFieldStats = p.Results.placeFieldStats;
UIDs = p.Results.UIDs;

% Get session info
basename = bz_BasenameFromBasepath(basepath);
load([basepath filesep basename '.sessionInfo.mat']);
% Default firingMapsAvg
if isempty(firingMaps)
    firingMaps = load([basepath filesep basename '.firingMapsAvg.cellinfo.mat']);
    firingMaps = firingMaps.firingMaps;
end
% Default firingMapsAvg
if isempty(placeFieldStats)
    placeFieldStats = load([basepath filesep basename '.placeFields.cellinfo.mat']);
    placeFieldStats = placeFieldStats.placeFieldStats;
end
saveMat = p.Results.saveMat;

%% Find place field template

% Number of Units
nUnits = length(placeFieldStats.mapStats);
% Number of conditions
nCond = length(placeFieldStats.mapStats{1});

% Default UIDs
if isempty(UIDs)
    UIDs = ones(nUnits);
end

% Template by center of biggest place field 
for c = 1:nCond
    for unit = 1:nUnits
        if UIDs(unit)
            templatePF{c}(unit,1) = placeFieldStats.mapStats{unit}{c}.x(1);
            templatePF{c}(unit,2) = placeFieldStats.UID(unit);
            templatePF{c}(unit,3) = unit;
        end
    end
end
% Sort units by place field, condition-independent
for c = 1:nCond
    templatePF{c} = sortrows(templatePF{c});
    templatePF{c}(isnan(templatePF{c}),:) = [];
end

% Template by center of Center of Mass
for c = 1:nCond
    for unit = 1:nUnits
        if UIDs(unit)
            x = 1:length( firingMaps.rateMaps{1}{1});
            rate = firingMaps.rateMaps{unit}{c};
            COM = sum(rate.*x) / sum(rate);
            templateCOM{c}(unit,1) = COM;
            templateCOM{c}(unit,2) = placeFieldStats.UID(unit);
            templateCOM{c}(unit,3) = unit;
        end
    end
end
% Sort units by Center of Mass, condition-independent
for c = 1:nCond
    templateCOM{c} = sortrows(templateCOM{c});
end

% Structure
placeFieldTemplate.Peak = templatePF;
placeFieldTemplate.CenterofMass = templateCOM;
% Save
if saveMat
   save([basepath filesep basename '.placeFieldTemplate.mat'], 'placeFieldTemplate'); 
end



end

