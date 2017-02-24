function [n,filters] = plx_chan_filters(filename)
% plx_chan_filters(filename): read channel filter settings for each spike channel from a .plx or .pl2 file
%
% [n,filters] = plx_chan_filters(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   filter - array of filter values (0 or 1)
%   n - number of channels

if nargin ~= 1
    error 'expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    n = numel(pl2.SpikeChannels);
    if n > 0
        filters = zeros(n,1);
    end
    return
end

[n,filters] = mexPlex(10, filename);

end