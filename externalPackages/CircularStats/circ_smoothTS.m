function [ts_smooth] = circ_smoothTS(varargin)
% USAGE
% [ts_smooth] = circ_smoothTS(varargin)
%
% Given a timeseries input of circular data, this function smooths over the
% desired number of bins and returns a circularly smooth timeseries
%
% INPUTS
%   
%   ts         a time series of circular data
%   nBins      integer number of bins you would like to smooth over
%   method     string argument that determines the smoothing method,
%              options are 'median', 'mean' [default: 'median']
%   exclude    vector of values in ts to exclude from smoothing (can be
%              used to exclude 0 values)
%
% OUTPUTS
%
%   ts_smooth  a time series vector of smoothed circular data
%
% HELP
%
% Written by David Tingley, 2017
% TODO error handling when nBins > 1/2 ts

p = inputParser;
addRequired(p,'ts',@isvector);
addRequired(p,'nBins',@isnumeric);
addParameter(p,'method','median',@isstr)
addParameter(p,'exclude',[],@isvector);

parse(p,varargin{:});

ts = p.Results.ts;
if size(ts,1) == 1;
    ts = ts';
end

nBins = p.Results.nBins;
method = p.Results.method;
exclude = p.Results.exclude;

if length(exclude) == length(ts)
   ts_smooth = ts;
   return
end
if nBins == 1
    ts_smooth = ts;
    return
end


if ~isempty(exclude)
    list = find(ts==exclude);
    ts(list) = nan;
end

exclude = find(isnan(ts));
for i = 1:length(ts)
    if i <= nBins
        ind = (1:i+floor(nBins/2)); 
    elseif i >= length(ts) - nBins
        ind = (i-floor(nBins/2):length(ts));
    else
        ind = (i-floor(nBins/2):i+floor(nBins/2));
    end
    
    [loc] = ~ismember(ind,exclude);
    if ~isempty(loc) & sum(loc) > 0
        switch method
            case 'median'
                ts_smooth(i) = circ_median(ts(ind(loc)));
            case 'mean'
                ts_smooth(i) = circ_mean(ts(ind(loc)));
            case 'guassian'

            case 'interp' 
                
        end
    else
        ts_smooth(i) = ts(i);    
    end
    
end

ts_smooth(isnan(ts_smooth)) = 0;








