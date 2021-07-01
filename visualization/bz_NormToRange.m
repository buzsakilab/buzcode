function [ normdata ] = bz_NormToRange(data,range,databounds)
%[ normdata ] = bz_NormToRange(data,range,databounds) normalizes some data to
%fit in some range.
%
%INPUTS
%   data        the data you want to normalize
%   range       [min max] you would like to normalize it to.
%               use 'ylim' to normalize to min/max of current plot (default)
%   databounds  [min max] of the data (optional)
%
%OUTPUTS
%   normdata    normalized data
%
%DLevenstein 2019
%%
if isempty(data); normdata = []; return; end

if ~exist('range','var') || strcmp(range,'ylim')
    range = get(gca,'ylim');
elseif length(range)==1
    bounds = get(gca,'ylim');
    bounds(2) = bounds(1)+range.*diff(bounds);
    range=bounds;
end


if ~exist('databounds','var')
    databounds(1) = min(data(:)); databounds(2) = max(data(:));
end


dataspan = diff(databounds);
rangespan = diff(range);

normdata = (data-databounds(1))./dataspan;
normdata = normdata.*rangespan+range(1);

end

