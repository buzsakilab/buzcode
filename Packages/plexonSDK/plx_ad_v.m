function [adfreq, n, ts, fn, ad] = plx_ad_v(filename, channel)
% plx_ad_v(filename, channel): read all a/d data for the specified channel from a .plx or .pl2 file
%              returns a/d values in millivolts
%
% [adfreq, n, ts, fn, ad] = plx_ad_v(filename, channel)
% [adfreq, n, ts, fn, ad] = plx_ad_v(filename, 0)
% [adfreq, n, ts, fn, ad] = plx_ad_v(filename, 'FP01')
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - 0-based channel number or channel name
%
%           a/d data come in fragments. Each fragment has a timestamp
%           and a number of a/d data points. The timestamp corresponds to
%           the time of recording of the first a/d value in this fragment.
%           All the data values stored in the vector ad.
% 
% OUTPUT:
%   adfreq - digitization frequency for this channel
%   n - total number of data points 
%   ts - array of fragment timestamps (one timestamp per fragment, in seconds)
%   fn - number of data points in each fragment
%   ad - array of a/d values converted to millivolts

adfreq = 0;
n = 0;
ts = -1;
fn = -1;
ad = -1;

if nargin ~= 2
    error 'Expected 2 input arguments';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    % return data from .pl2 file
    pl2Channel = channel;
    % PL2 uses 1-based channel numbers
    if ischar(channel) == 0
        pl2Channel = channel+1;
    end
    pl2ad = PL2Ad(filename, pl2Channel);
    if numel(pl2ad.Values) > 0
        adfreq = pl2ad.ADFreq;
        n = numel(pl2ad.Values);
        ts = pl2ad.FragTs;
        fn = pl2ad.FragCounts;
        ad = pl2ad.Values;
    end
    return;
end

channelNumber = plx_ad_resolve_channel(filename, channel);
if channelNumber == -1
    fprintf('\n plx_ad_v: no header for the specified A/D channel.');
    fprintf('\n    use plx_ad_info(filename) to print the list of a/d channels\n');
    return
end

[adfreq, n, ts, fn, ad] = mexPlex(17, filename, channelNumber);

end