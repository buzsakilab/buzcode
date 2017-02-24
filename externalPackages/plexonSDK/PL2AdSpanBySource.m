function [ad] = PL2AdSpanBySource(filename, source, channelInSource, startCount, endCount)
% PL2AdSpanBySource(filename, source, channelInSource, startCount, endCount): read a span of a/d data from a .pl2 file
%                returns a/d values in millivolts
%
% examples:
%
% ad = PL2AdSpanBySource(filename, source, channelInSource, startCount, endCount)
% ad = PL2AdSpanBySource(filename, 'FP', 1, 10000)
% ad = PL2AdSpanBySource(filename, 8, 1, 1000, 2000)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   source - source name or source numeric id
%   channel -  1-based channel number within source
%           To print the list of analog channels, use these commands:
%
%           pl2 = PL2GetFileIndex(filename);
%           PL2Print(pl2.AnalogChannels);
%
%   startCount - index of first sample to fetch
%   endCount - index of last sample to fetch
%
% OUTPUT:
%   ad.Values - array of a/d values in millivolts
%   ad.ADFreq - digitization frequency for this channel (samples per second)

ad.Values = [];
ad.ADFreq = 0;

if nargin ~= 5
    error 'expected 5 input arguments';
end

filename = internalPL2ResolveFilename(filename);
pl2 = PL2GetFileIndex(filename);
channelNumber = internalPL2ResolveChannelBySource(pl2.AnalogChannels, source, channelInSource);
ad = Pl2AdSpan(filename, channelNumber, startCount, endCount);

end