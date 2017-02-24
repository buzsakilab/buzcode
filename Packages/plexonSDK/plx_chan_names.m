function [n,names] = plx_chan_names(filename)
% plx_chan_names(filename): read name for each spike channel from a .plx or .pl2 file
%
% [n,names] = plx_chan_names(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   names - array of channel name strings
%   n - number of channels

if nargin ~= 1
    error 'expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    n = numel(pl2.SpikeChannels);
    if n > 0
        for i=1:n
            if i == 1
                names = pl2.SpikeChannels{i}.Name;
            else
                names = char(names, pl2.SpikeChannels{i}.Name);
            end
        end
        names = internalPL2TrimNames(names);
    end
    return;
end

[n,names] = mexPlex(14, filename);

end