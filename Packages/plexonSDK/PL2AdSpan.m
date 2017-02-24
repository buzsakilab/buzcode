function [ad] = PL2AdSpan(filename, channel, startCount, endCount)
% PL2AdSpan(filename, channel, startCount, endCount): read a span of a/d data from a .pl2 file
%                returns a/d values in millivolts
%
% examples:
%
% ad = PL2AdSpan(filename, channel, startCount, endCount)
% ad = PL2AdSpan(filename, 48, 1, 10000)
% ad = PL2AdSpan(filename, 'FP02', 1000, 2000)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel -  channel name or 1-based channel number
%           1-based channel number is the index of the a/d channel in pl2.AnalogChannels cell array
%           shown in the Name column of the output produced by PL2Print(pl2.AnalogChannels)
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

[ad] = internalPL2AdSpan(filename, channel, startCount, endCount);

end