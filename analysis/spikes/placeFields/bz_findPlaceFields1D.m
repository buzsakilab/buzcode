function [mapStats] = bz_findPlaceFields1D(firingMaps)
%   [stats] = bz_findPlaceFields1D(firingMaps)
%   Find place fields from 1D firing maps 

%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'       number of bins (default = 50)
%     'minTime'     minimum time spent in each bin (in s, default = 0)
%     'maxGap'      z values recorded during time gaps between successive (x,y)
%                   samples exceeding this threshold (e.g. undetects) will not
%                   be interpolated; also, such long gaps in (x,y) sampling
%                   will be clipped to 'maxGap' to compute the occupancy map
%                   (default = 0.100 s)
%     'type'        'linear' for linear data, 'circular' for angular data
%                   (default 'linear')
%     'threshold'   values above threshold*peak belong to the field
%                   (default = 0.2)
%     'minSize'     fields smaller than this size are considered spurious
%                   and ignored (default = 10)
%     'minPeak'     peaks smaller than this size are considered spurious
%                   and ignored (default = 1)
%     'verbose'     display processing information (default = 'off')
%    =========================================================================

%% Parse arguments


%% Calculate
for unit = 1:length(spikes.times)
    for c = 1:conditions
        mapStats{unit}{c} = MapStats(firingMaps.rateMaps{unit}{c});
    end
end


%% Write output



end

