function [n, samplecounts] = plx_adchan_samplecounts(filename)
% plx_adchan_samplecounts(filename): read the per-channel sample counts for analog channels from a .plx or .pl2 file
%
% [n, samplecounts] = plx_adchan_samplecounts(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   n - number of channels
%   samplecounts - array of sample counts

if nargin ~= 1
    error 'expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    n = numel(pl2.AnalogChannels);
    samplecounts = zeros(n,1);
    for i=1:n
        samplecounts(i,1) = pl2.AnalogChannels{i}.NumValues;
    end
    return;
end

[n,samplecounts] = mexPlex(23, filename);

end