function [blockInfo] = PL2ReadFirstDataBlock(filename)
% PL2ReadFirstDataBlock(filename): read first data block, return block info
%
% example:
%
% blockInfo = PL2ReadFirstDataBlock(filename)
% while blockInfo.Type > 0
%    % process data from the block
%    % ...
%    % read the next data block
%    blockInfo = internalPL2ReadNextDataBlock(blockInfo);
% end
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   blockInfo - information about the first data block
%
%   blockInfo.Type - type of data block
%                    0 - not a data block, did not find data in the file
%                    1 - spike data type
%                    2 - analog data type
%                    3 - digital event data type
%                    4 - start/stop event data type
%   blockInfo.Source - block source
%   blockInfo.Channel - block channel inside source
%   if blockInfo.Type is 1 (spike), SpikeData structure is filled
%       blockInfo.SpikeData.WaveLength - number of samples in each waveform
%       blockInfo.SpikeData.NumSpikes - number of spikes 
%       blockInfo.SpikeData.Timestamps  - spike timestamps in seconds
%       blockInfo.SpikeData.Units - spike units 
%       blockInfo.SpikeData.Waveforms - spike waveforms in milliVolts
%   if blockInfo.Type is 2 (analog), AnalogData structure is filled
%       blockInfo.AnalogData.NumSamples - number of analog samples
%       blockInfo.AnalogData.Timestamp - timestamps of the first sample (in seconds)
%       blockInfo.AnalogData.Values - sample values in milliVolts
%   if blockInfo.Type is 3 (event), EventData structure is filled
%       blockInfo.EventData.NumEvents - number of events
%       blockInfo.EventData.Timestamps - event timestamps in seconds
%       blockInfo.EventData.Values - event values
%   if blockInfo.Type is 4 (start/stop event), StartStopEventData structure is filled
%       blockInfo.StartStopEventData.NumEvents - number of events
%       blockInfo.StartStopEventData.Timestamps - event timestamps in seconds
%       blockInfo.StartStopEventData.Values - event values

blockInfo = internalPL2ReadFirstDataBlock(filename);

end