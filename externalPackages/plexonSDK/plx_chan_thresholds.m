function [n,thresholds] = plx_chan_thresholds(filename)
% plx_chan_thresholds(filename): read channel thresholds from a .plx or .pl2 file
%
% [n,thresholds] = plx_chan_thresholds(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   thresholds - array of tresholds, expressed in raw A/D counts
%   n - number of channel

if nargin ~= 1
    error 'expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    n = numel(pl2.SpikeChannels);
    if n > 0
        thresholds = zeros(n,1);
        for i=1:n
            thresholds(i,1) = pl2.SpikeChannels{i}.Threshold;
        end
    end
    return;
end

[n,thresholds] = mexPlex(9, filename);

end