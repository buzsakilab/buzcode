function [ad] = PL2Ad(filename, channel)
% PL2Ad(filename, channel): read all a/d data for the specified channel from a .pl2 file
%                returns a/d values in millivolts
%
% examples:
%
% ad = PL2Ad(filename, channel)
% ad = PL2Ad(filename, 2)
% ad = PL2Ad(filename, 'FP02')
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

[ad] = internalPL2Ad(filename, channel);

end