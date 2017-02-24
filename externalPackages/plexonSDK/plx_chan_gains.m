function [n,gains] = plx_chan_gains(filename)
% plx_chan_gains(filename): read channel gains from .plx or .pl2 file
%
% [gains] = plx_chan_gains(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%  gains - array of total gains
%   n - number of channels

if nargin ~= 1
    error 'expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    n = numel(pl2.SpikeChannels);
    if n > 0
        gains = zeros(n,1);
        for i=1:n
            gains(i,1) = pl2.SpikeChannels{i}.TotalGain;
        end
    end
    return;
end

[n,gains] = mexPlex(8, filename);

end