function [] = plx_ad_info(filename)
% plx_ad_info(filename): prints .plx file a/d channel info
%
% plx_ad_info(filename)
%
% INPUT:
%   filename - .plx file name. Will use file open dialog if filename is empty string
%
% OUTPUT:
%   (none)

if nargin ~= 1
    error 'expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    PL2Print(pl2.AnalogChannels);
    return
end

[n, names] = plx_adchan_names(filename);
[n, raw] = plx_ad_chanmap(filename);
[n, samplecounts] = plx_adchan_samplecounts(filename);
[n, freqs] = plx_adchan_freqs(filename);
fprintf(' %18.18s', 'RawChannelNum');
fprintf(' %18.18s', 'Name');
fprintf(' %18.18s', 'SampleCount');
fprintf(' %18.18s', 'SamplesPerSecond');
fprintf('\n');
    
for i=1:n
    fprintf(' %18d', raw(i));
    fprintf(' %18.18s', names(i,:));
    fprintf(' %18d', samplecounts(i));
    fprintf(' %18.10g', freqs(i));
    fprintf('\n');
end

end