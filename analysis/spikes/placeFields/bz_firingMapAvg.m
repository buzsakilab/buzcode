function [firingMaps] = bz_firingMapAvg(positions,spikes,varargin)

% USAGE
% [firingMaps] = bz_firingMapAvg(positions,spikes,varargin)
% Calculates averaged firing maps for each cell in 1D or 2D enviroments
%
% INPUTS
%
%   spikes    - buzcode format .cellinfo. struct with the following fields
%               .times 
%   positions - [t x y ] or [t x] position matrix or
%               cell with several of these matrices (for different conditions)
%      or
%   behavior  - buzcode format behavior struct - NOT YET IMPLEMENTED
%   <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'       number of bins (default = 50)
%     'minTime'     minimum time spent in each bin (in s, default = 0)
%     'mode'        'interpolate' to interpolate missing points (< minTime),
%                   or 'discard' to discard them (default)
%     'maxDistance' maximal distance for interpolation (default = 5)
%     'maxGap'      z values recorded during time gaps between successive (x,y)
%                   samples exceeding this threshold (e.g. undetects) will not
%                   be interpolated; also, such long gaps in (x,y) sampling
%                   will be clipped to 'maxGap' to compute the occupancy map
%                   (default = 0.100 s)
%     'type'        'linear' for linear data, 'circular' for angular data
%                   (default 'linear')
%      saveMat   - logical (default: false) that saves firingMaps file
%      CellInspector   - logical (default: false) that creates an otuput
%                   compatible with CellInspector

%
%
% OUTPUT
%
%   firingMaps - cellinfo struct with the following fields
%                .rateMaps              gaussian filtered rates
%                .rateMaps_unsmooth     raw rate data
%                .rateMaps_box          box filtered rates
%                .countMaps             raw spike count data
%                .occuMaps              position occupancy data
%
% Antonio FR, 10/2019

%% parse inputs
p=inputParser;
addParameter(p,'smooth',2,@isnumeric);
addParameter(p,'nBins',50,@isnumeric);
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'type','linear',@isstr);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'CellInspector',false,@islogical);
addParameter(p,'mode','discard',@isstr);
addParameter(p,'maxDistance',5,@isnumeric);

parse(p,varargin{:});
smooth = p.Results.smooth;
nBins = p.Results.nBins;
maxGap = p.Results.maxGap;
minTime = p.Results.minTime;
type = p.Results.type;
saveMat = p.Results.saveMat;
CellInspector = p.Results.CellInspector;
mode = p.Results.mode;
maxDistance = p.Results.maxDistance;

% number of conditions
  if iscell(positions)
     conditions = length(positions); 
  elseif isvector(positions)
     conditions = 1;
  end
  %%% TODO: conditions label
  
%% Calculate

% speed threshold 



% get firign rate maps
for unit = 1:length(spikes.times)
    for c = 1:conditions
        map{unit}{c} = Map(positions{c},spikes.times{unit},'smooth',smooth,'minTime',minTime,...
            'nBins',nBins,'maxGap',maxGap,'mode',mode,'type',type,'maxDistance',maxDistance);
    end
end

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
firingMaps.UID = spikes.UID;
firingMaps.sessionName = spikes.sessionName;
try
firingMaps.region = spikes.region; 
catch
   %warning('spikes.region is missing') 
end

firingMaps.params.smooth = smooth;
firingMaps.params.nBins = nBins;
firingMaps.params.maxGap = maxGap;

for unit = 1:length(spikes.times)
    for c = 1:conditions
    firingMaps.rateMaps{unit,1}{c} = map{unit}{c}.z;
    firingMaps.countMaps{unit,1}{c} = map{unit}{c}.count;
    firingMaps.occupancy{unit,1}{c} = map{unit}{c}.time;
    end
end

if saveMat
   save([firingMaps.sessionName '.firingMapsAvg.cellinfo.mat'],'firingMaps'); 
end

% output for cell inspector - only works for 1D
if CellInspector
    for unit = 1:length(spikes.times)
        for c = 1:conditions
            firingRateMap.map{unit}(:,c) = map{unit}{c}.z;
        end
    end
    firingRateMap.x_bins = 1:1:lenght(map{1}{1}.z);
    firingRateMap.state_labels = nan(1,conditions);
            
    save([firingMaps.sessionName '.firingRateMap.mat'],'firingMaps'); 

end

%% Add option to plot

   
end
