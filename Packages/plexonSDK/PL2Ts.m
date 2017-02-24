function [ts] = PL2Ts( filename, channel, unit )
% PL2Ts( filename, channel, unit ): read spike timestamps from a .pl2 file
%
% ts = PL2Ts( filename, channel, unit );
% ts = PL2Ts( filename, 'SPK01', 1 );
% ts = PL2Ts( filename, 2, 3 );
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - channel name or 1-based channel number
%           1-based channel number is the index of the spike channel in pl2.SpikeChannels cell array
%           shown in the Name column of the output produced by PL2Print(pl2.SpikeChannels)
%           To print the list of spike channels, use these commands:
%
%           pl2 = PL2GetFileIndex(filename);
%           PL2Print(pl2.SpikeChannels);
%
%   unit  - unit number (0 - unsorted, 1-255 sorted units)
%
% OUTPUT:
%   ts - array of timestamps (in seconds)

ts = internalPL2Ts( filename, channel, unit );

end