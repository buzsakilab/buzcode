function [wave] = PL2Waves(filename, channel, unit)
% PL2Waves( filename, channel, unit ): read spike timestamps and waveforms from a .pl2 file
%                returns waveform values in millivolts
%
% wave = PL2Waves( filename, channel, unit )
% wave = PL2Waves( filename, 'SPK01', 1 );
% wave = PL2Waves( filename, 2, 3 );
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
%   unit  - unit number (0 - unsorted, 1-254 sorted units, 255 - invalidated spikes)
%
% OUTPUT:
%   wave.NumPointsWave - number of points in each waveform
%   wave.Ts - array of timestamps (in seconds) 
%   wave.Waves - matrix of waveforms, each row is a waveform in millivolts

wave = internalPL2Waves(filename, channel, unit);

end