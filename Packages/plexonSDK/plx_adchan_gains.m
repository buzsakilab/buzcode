function [n,gains] = plx_adchan_gains(filename)
% plx_adchan_gains(filename): read analog channel gains from .plx or .pl2 file
%
% [n,gains] = plx_adchan_gains(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%  gains - array of total gains
%  n - number of channels

if nargin ~= 1
    error 'expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    n = numel(pl2.AnalogChannels);
    gains = zeros(n,1);
    for i=1:n
        gains(i,1) = pl2.AnalogChannels{i}.TotalGain;
    end
    return;
end

[n,gains] = mexPlex(11, filename);

end