function [n, freqs] = plx_adchan_freqs(filename)
% plx_adchan_freq(filename): read the per-channel frequencies for analog channels from a .plx or .pl2 file
%
% [n, freqs] = plx_adchan_freq(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   freqs - array of frequencies
%   n - number of channels

if nargin ~= 1
    error 'expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    n = numel(pl2.AnalogChannels);
    if n > 0
        freqs = zeros(n,1);
        for i=1:n
            freqs(i,1) = pl2.AnalogChannels{i}.SamplesPerSecond;
        end
    end
    return;
end

[n,freqs] = mexPlex(12, filename);

end