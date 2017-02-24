function [adfreq, n, ad] = plx_ad_span(filename, channel, startCount, endCount)
% plx_ad_span(filename, channel): read a span of a/d data from a .plx or .pl2 file
%
% [adfreq, n, ad] = plx_ad_span(filename, channel, startCount, endCount)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   startCount - index of first sample to fetch
%   endCount - index of last sample to fetch
%   channel - 0 - based channel number or channel name
%
% OUTPUT:
%   adfreq - digitization frequency for this channel
%   n - total number of data points 
%   ad - array of raw a/d values

adfreq = 0;
n = 0;
ad = -1;

if nargin ~= 4
    error 'expected 4 input arguments';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    % return data from .pl2 file
    pl2Channel = channel;
    % PL2 uses 1-based channel numbers
    if ischar(channel) == 0
        pl2Channel = channel+1;
    end
    pl2ad = PL2AdSpan(filename, pl2Channel, startCount, endCount);
    if numel(pl2ad.Values) > 0
        adfreq = pl2ad.ADFreq;
        n = numel(pl2ad.Values);
        ad = pl2ad.Values;
        % convert millivolts to raw a/d values
        pl2 = PL2GetFileIndex(filename);
        channelNumber = internalPL2ResolveChannel(pl2.AnalogChannels, pl2Channel);
        ad = round(ad / (pl2.AnalogChannels{channelNumber}.CoeffToConvertToUnits * 1000));  
    end
    return
end

% return plx data
channelNumber = plx_ad_resolve_channel(filename, channel);
if channelNumber == -1
    fprintf('\n plx_ad_span: no header for the specified A/D channel.');
    fprintf('\n    use plx_ad_info(filename) to print the list of a/d channels\n');
    return
end

[adfreq, n, ad] = mexPlex(7, filename, channelNumber, startCount, endCount);
end