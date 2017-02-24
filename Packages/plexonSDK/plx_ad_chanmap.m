function  [n, adchans] = plx_ad_chanmap(filename)
% plx_ad_chanmap(filename): return map of raw continuous channel numbers for each channel
%            for the specified .plx or .pl2 file
%
% [n, adchans] = plx_ad_chanmap(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   n - number of continuous channels
%   adchans - 1 x n array of continuous channel numbers
%
% Normally, there is one channel entry in the .plx for for each raw continuous channel,
% so the mapping is trivial adchans[i] = i-1 (because continuous channels start at 0).
% However, for certain .plx files saved in some ways from OFS (notably after
% loading data files from other vendors), the mapping can be more complex.
% E.g. there may be only 2 non-empty channels in a .plx file, but those channels
% correspond to raw channel numbers 7 and 34. So in this case NChans = 2, 
% and adchans[1] = 7, adchans[2] = 34.
% The plx_ routines that return arrays always return arrays of size NChans. However,
% routines that take channels numbers as arguments always expect the raw  
% channel number.  So in the above example, to get the data from  
% the second channel, use
%   [adfreq, n, ts, fn, ad] = plx_ad(filename, adchans[2])

n = 0;
adchans = -1;

if nargin ~= 1
    error 'Expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    n = numel(pl2.AnalogChannels);
    adchans = 0:(n-1);
    return;
end

[n, adchans] = mexPlex(27,filename);

end