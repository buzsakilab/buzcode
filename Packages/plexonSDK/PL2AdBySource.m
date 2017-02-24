function [ad] = PL2AdBySource(filename, source, channelInSource)
% PL2AdBySource(filename, source, channelInSource): read all a/d data for the specified channel from a .pl2 file
%                returns a/d values in millivolts
%
% examples:
%
% ad = PL2AdBySource(filename, source, channelInSource)
% ad = PL2AdBySource(filename, 'FP', 1)
% ad = PL2AdBySource(filename, 8, 3)
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
% OUTPUT:
%           a/d data come in fragments. Each fragment has a timestamp
%           and a number of a/d data points. The timestamp corresponds to
%           the time of recording of the first a/d value in this fragment.
%           All the data values are stored in the vector ad.Values.
%
%   ad.Values - array of a/d values in milliVolts
%   ad.FragTs - array of fragment timestamps (one timestamp per fragment, in seconds)
%   ad.FragCounts - number of data points in each fragment
%   ad.ADFreq - digitization frequency for this channel (samples per second)

ad.Values = [];
ad.FragTs = [];
ad.FragCounts = [];
ad.ADFreq = 0;

if nargin ~= 3
    error 'expected 3 input arguments';
end

filename = internalPL2ResolveFilename(filename);
pl2 = PL2GetFileIndex(filename);
channelNumber = internalPL2ResolveChannelBySource(pl2.AnalogChannels, source, channelInSource);
ad = PL2Ad(filename, channelNumber);

end