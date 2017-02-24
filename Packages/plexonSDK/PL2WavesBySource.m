function [wave] = PL2WavesBySource(filename, source, channelInSource, unit)
% PL2WavesBySource( filename, source, channelInSource, unit ): read spike timestamps and waveforms from a .pl2 file
%                returns waveform values in millivolts
%
% wave = PL2WavesBySource( filename, source, channelInSource, unit )
% wave = PL2WavesBySource( filename, 'SPK', 1, 1 );
% wave = PL2WavesBySource( filename, 7, 2, 3 );
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
%   wave.NumPointsWave - number of points in each waveform
%   wave.Ts - array of timestamps (in seconds) 
%   wave.Waves - matrix of waveforms, each row is a waveform in millivolts

wave.NumPointsWave = 0;
wave.Ts = [];
wave.Waves = [];

if nargin ~= 4
    error 'expected 4 input arguments';
end

filename = internalPL2ResolveFilename(filename);
pl2 = PL2GetFileIndex(filename);
channelNumber = internalPL2ResolveChannelBySource(pl2.SpikeChannels, source, channelInSource);
wave = PL2Waves(filename, channelNumber, unit);

end