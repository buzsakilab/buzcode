function [ad] = PL2AdTimeSpanBySource(filename, source, channelInSource, startTime, endTime)
% PL2AdTimeSpanBySource(filename, source, channelInSource, startTime, endTime): read a span of a/d data from a .pl2 file
%                returns a/d values in millivolts
%
% examples:
%
% ad = PL2AdTimeSpanBySource(filename, source, channelInSource, startTime, endTime)
% ad = PL2AdTimeSpanBySource(filename, 'FP', 2, 0.5, 60)
% ad = PL2AdTimeSpanBySource(filename, 8, 1, 0, 10.5)
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
%   startTime - start time in seconds
%   endTime - end time in seconds
%
% OUTPUT:
%           a/d data come in fragments. Each fragment has a timestamp
%           and a number of a/d data points. The timestamp corresponds to
%           the time of recording of the first a/d value in this fragment.
%           All the data values are stored in the vector ad.Values.
%
%   ad.Values - array of a/d values in millivolts
%   ad.FragTs - array of fragment timestamps (one timestamp per fragment, in seconds)
%   ad.FragCounts - number of data points in each fragment
%   ad.ADFreq - digitization frequency for this channel (samples per second)

ad.Values = [];
ad.FragTs = [];
ad.FragCounts = [];
ad.ADFreq = 0;

if nargin ~= 5
    error 'expected 5 input arguments';
end

filename = internalPL2ResolveFilename(filename);
pl2 = PL2GetFileIndex(filename);
channelNumber = internalPL2ResolveChannelBySource(pl2.AnalogChannels, source, channelInSource);
ad = PL2AdTimeSpan(filename, channelNumber, startTime, endTime);

end