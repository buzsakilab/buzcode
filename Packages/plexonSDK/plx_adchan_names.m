function [n,names] = plx_adchan_names(filename)
% plx_adchan_names(filename): gets the names of a/d channels in a .plx or .pl2 file
%
% [n,names] = plx_adchan_names(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   names - array of a/d channel name strings
%   n - number of channels

if nargin ~= 1
    error 'expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    n = numel(pl2.AnalogChannels);
    for i=1:numel(pl2.AnalogChannels)
        if i == 1
            names = pl2.AnalogChannels{i}.Name;
        else
            names = char(names, pl2.AnalogChannels{i}.Name);
        end
    end
    names = internalPL2TrimNames(names);
    return;
end

[n,names] = mexPlex(15, filename);

end