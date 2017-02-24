function  [n,dspchans] = plx_chanmap(filename)
% plx_chanmap(filename): return map of raw DSP channel numbers for each channel
%
% [n,dspchans] = plx_chanmap(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   n - number of spike channels
%   dspchans - 1 x n array of DSP channel numbers
%
% Normally, there is one channel entry in the .plx for for each raw DSP channel,
% so the mapping is trivial dspchans[i] = i.
% However, for certain .plx files saved in some ways from OFS (notably after
% loading data files from other vendors), the mapping can be more complex.
% E.g. there may be only 2 non-empty channels in a .plx file, but those channels
% correspond to raw DSP channel numbers 7 and 34. So in this case NChans = 2, 
% and dspchans[1] = 7, dspchans[2] = 34.
% The plx_ routines that return arrays always return arrays of size NChans. However,
% routines that take channels numbers as arguments always expect the raw DSP 
% channel number.  So in the above example, to get the timestamps from unit 4 on 
% the second channel, use
%   [n,ts] = plx_ts(filename, dspchans[2], 4 );

if nargin ~= 1
    error 'expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    n = numel(pl2.SpikeChannels);
    dspchans = 1:n;
    return;
end

[n,dspchans] = mexPlex(26, filename);

end