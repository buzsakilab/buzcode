function [ts] = PL2TsBySource( filename, source, channelInSource, unit )
% PL2TsBySource( filename, source, channelInSource, unit ): read spike timestamps from a .pl2 file
%
% ts = PL2TsBySource( filename, source, channelInSource, unit );
% ts = PL2TsBySource( filename, 'SPK', 1, 2 );
% ts = PL2TsBySource( filename, 7, 1, 3 );
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   source - source name or source numeric id
%   channel -  1-based channel number within source
%           To print the list of spike channels, use these commands:
%
%           pl2 = PL2GetFileIndex(filename);
%           PL2Print(pl2.SpikeChannels);
%
%   unit  - unit number (0 - unsorted, 1-255 sorted units)
%
% OUTPUT:
%   ts - array of timestamps (in seconds)

n = 0;
ts = [];

if nargin ~= 4
    error 'expected 4 input arguments';
end

filename = internalPL2ResolveFilename(filename);
pl2 = PL2GetFileIndex(filename);
channelNumber = internalPL2ResolveChannelBySource(pl2.SpikeChannels, source, channelInSource);
ts = Pl2Ts(filename, channelNumber, unit);

end