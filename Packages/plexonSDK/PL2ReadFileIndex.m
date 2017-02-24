function [ pl2 ] = PL2ReadFileIndex( filename )
% PL2ReadFileIndex( filename ): reads pl2 file index. 
%      Displays error if the file does not exist.
%      This function is not called directly by PL2* functions. 
%      Use PL2GetFileIndex(filename) instead. 
%      PL2GetFileIndex caches file indexes so it works faster if you run multiple PL* functions on the same file.                        
%
% [pl2] = PL2ReadFileIndex( filename )
%
% INPUT:
%   filename - pl2 file path. 
%
% OUTPUT:
%   pl2 - pl2 file index object 
%       
%   pl2.FilePath - path of the file
%   pl2.FileLength - file length in bytes
%
%   pl2.AnalogChannels - array of analog channel structures
%   pl2.AnalogChannels{i}.PlxChannel - 0-based .plx channel number 
%   pl2.AnalogChannels{i}.NumValues - number of analog values
%   pl2.AnalogChannels{i}.Name - channel name
%   pl2.AnalogChannels{i}.SourceName - source name
%   pl2.AnalogChannels{i}.Source - source numeric ID
%   pl2.AnalogChannels{i}.Channel - 1-based channel number within source 
%   pl2.AnalogChannels{i}.Enabled - (zero or one) flag indicating if channel is enabled 
%   pl2.AnalogChannels{i}.RecordingEnabled - (zero or one) flag indicating if channel recording is enabled 
%   pl2.AnalogChannels{i}.Units - units for channel values  
%   pl2.AnalogChannels{i}.SamplesPerSecond - channel digitizing rate
%   pl2.AnalogChannels{i}.CoeffToConvertToUnits - coefficient to convert raw a/d values to units
%   pl2.AnalogChannels{i}.SourceTrodality - source trodality
%   pl2.AnalogChannels{i}.Trode - 1-based trode within source
%   pl2.AnalogChannels{i}.ChannelInTrode - 1-based channel within trode
%   pl2.AnalogChannels{i}.NumberOfChannelsInSource - number of channels in source
%   pl2.AnalogChannels{i}.DeviceId - source device id
%   pl2.AnalogChannels{i}.NumberOfChannelsInDevice - number of channels in device
%   pl2.AnalogChannels{i}.SourceDeviceName - source device name
%   pl2.AnalogChannels{i}.ProbeDeviceName - device name of the probe that is used to record the channel
%   pl2.AnalogChannels{i}.ProbeSourceId - source numeric ID of the probe
%   pl2.AnalogChannels{i}.ProbeSourceChannel - 1-based channel number of the probe within probe source
%   pl2.AnalogChannels{i}.ProbeDeviceId - device ID of the probe source
%   pl2.AnalogChannels{i}.ProbeDeviceChannel - 1-based channel number of the probe within probe device
%   pl2.AnalogChannels{i}.InputVoltageMinimum - a/d converter input voltage minimum in volts
%   pl2.AnalogChannels{i}.InputVoltageMaximum - a/d converter input voltage maximum in volts
%   pl2.AnalogChannels{i}.TotalGain - total signal gain before siginal is digitized in a/d converter
%
%   pl2.SpikeChannels - array of spike channel structures
%   pl2.SpikeChannels{i}.Name - channel name
%   pl2.SpikeChannels{i}.SourceName - source name
%   pl2.SpikeChannels{i}.Source - source numeric ID
%   pl2.SpikeChannels{i}.Channel - 1-based channel number within source 
%   pl2.SpikeChannels{i}.Enabled - (zero or one) flag indicating if channel is enabled 
%   pl2.SpikeChannels{i}.RecordingEnabled - (zero or one) flag indicating if channel recording is enabled 
%   pl2.SpikeChannels{i}.Units - units for channel values  
%   pl2.SpikeChannels{i}.SamplesPerSecond - channel digitizing rate
%   pl2.SpikeChannels{i}.CoeffToConvertToUnits - coefficient to convert raw a/d values to units
%   pl2.SpikeChannels{i}.SamplesPerSpike - number of samples in each waveform
%   pl2.SpikeChannels{i}.Threshold - threshold in raw a/d values
%   pl2.SpikeChannels{i}.PreThresholdSamples - number of camples before signal crosses the threshold
%   pl2.SpikeChannels{i}.SortEnabled - sort enabled flag
%   pl2.SpikeChannels{i}.SortMethod - sort method
%   pl2.SpikeChannels{i}.NumberOfUnits - number of units
%   pl2.SpikeChannels{i}.SortRangeStart - sort range start (samples)
%   pl2.SpikeChannels{i}.SortRangeEnd - sort range end (samples)
%   pl2.SpikeChannels{i}.UnitCounts - unit counts. UnitCounts(1) - number of unsorted spikes
%                                                  UnitCounts(2) - number of unit 1 spikes, etc.
%												   UnitCounts(256) - number of invalidated spikes
%   pl2.SpikeChannels{i}.SourceTrodality - source trodality
%   pl2.SpikeChannels{i}.Trode - 1-based trode within source
%   pl2.SpikeChannels{i}.ChannelInTrode - 1-based channel within trode
%   pl2.SpikeChannels{i}.NumberOfChannelsInSource - number of channels in source
%   pl2.SpikeChannels{i}.DeviceId - source device id
%   pl2.SpikeChannels{i}.NumberOfChannelsInDevice - number of channels in device
%   pl2.SpikeChannels{i}.SourceDeviceName - source device name
%   pl2.SpikeChannels{i}.ProbeDeviceName - device name of the probe that is used to record the channel
%   pl2.SpikeChannels{i}.ProbeSourceId - source numeric ID of the probe
%   pl2.SpikeChannels{i}.ProbeSourceChannel - 1-based channel number of the probe within probe source
%   pl2.SpikeChannels{i}.ProbeDeviceId - device ID of the probe source
%   pl2.SpikeChannels{i}.ProbeDeviceChannel - 1-based channel number of the probe within probe device
%   pl2.SpikeChannels{i}.InputVoltageMinimum - a/d converter input voltage minimum in volts
%   pl2.SpikeChannels{i}.InputVoltageMaximum - a/d converter input voltage maximum in volts
%   pl2.SpikeChannels{i}.TotalGain - total signal gain before siginal is digitized in a/d converter
%   
%   pl2.EventChannels - array of event channel structures
%   pl2.EventChannels{i}.NumEvents - number of events
%   pl2.EventChannels{i}.Name - channel name
%   pl2.EventChannels{i}.SourceName - source name
%   pl2.EventChannels{i}.Source - source numeric ID
%   pl2.EventChannels{i}.Channel - 1-based channel number within source 
%   pl2.EventChannels{i}.Enabled - (zero or one) flag indicating if channel is enabled 
%   pl2.EventChannels{i}.RecordingEnabled - (zero or one) flag indicating if channel recording is enabled 
%   pl2.EventChannels{i}.SourceDeviceName - source device name
%   pl2.EventChannels{i}.NumberOfChannelsInSource - number of channels in source
%   pl2.EventChannels{i}.DeviceId - source device id
%   pl2.EventChannels{i}.NumberOfChannelsInDevice - number of channels in device
%   pl2.SpikeChannels{i}.SourceDeviceName - source device name
%   pl2.EventChannels{i}.PlxChannel - channel number in .plx file
%
%   pl2.StartStopChannel - start/stop channel structure
%   pl2.StartStopChannel.NumEvents - number of events
%
%   pl2.Version - pl2 file version
%   pl2.Version.Major - major file version
%   pl2.Version.Minor - minor file version
%   pl2.Version.Bugfix - bugfix file version
%
%   pl2.StartRecordingTimeTicks - OmniPlex internal time of start of recording (in ticks)
%   pl2.DurationOfRecordingTicks - duration of recording (in ticks)
%   pl2.CreatorComment - file comment
%   pl2.CreatorSoftwareName - name of the program that created data file
%   pl2.CreatorSoftwareVersion - version of the program that created data file
%
%   pl2.CreatorDateTime - date and time when the file was created
%   pl2.CreatorDateTime.Second - seconds after the minute - [0,59]
%   pl2.CreatorDateTime.Minute - minutes after the hour - [0,59]
%   pl2.CreatorDateTime.Hour - hours since midnight - [0,23] 
%   pl2.CreatorDateTime.MonthDay - day of the month - [1,31]
%   pl2.CreatorDateTime.Month - months since January - [0,11]
%   pl2.CreatorDateTime.Year - years since 1900
%   pl2.CreatorDateTime.WeekDay - days since Sunday - [0,6]
%   pl2.CreatorDateTime.YearDay - days since January 1 - [0,365]
%   pl2.CreatorDateTime.IsDst - daylight savings time flag
%   pl2.CreatorDateTime.MilliSec - milliseconds after the second [0,999]
%
%   pl2.TimestampFrequency - file timestamp frequency in hertz
%   pl2.DurationOfRecordingSec - duration of recording (in seconds)
%
%   pl2.TotalNumberOfSpikeChannels - total number of spike channels
%   pl2.NumberOfRecordedSpikeChannels - number of recorded spike channels
%   pl2.TotalNumberOfAnalogChannels - total number of analog channels
%   pl2.NumberOfRecordedAnalogChannels - number of recorded analog channels
%   pl2.NumberOfEventChannels - number of event channels
%
%   pl2.MinimumTrodality - minimum trodality (number of channels in trode)
%   pl2.MaximumTrodality - maximum trodality (number of channels in trode)
%
%   pl2.NumberOfNonOmniPlexSources - number of non-OmniPlex sources (arbitrary binary data, can be used for CinePlex data)
%
%   pl2.ReprocessorComment - reprocessor comment
%   pl2.ReprocessorSoftwareName - reprocessor software name
%   pl2r.ReprocessorDateTime - reprocessor date and time
%
%   pl2.EventChanmap - .plx event channel map
%   pl2.EventCounts - .plx event counts

pl2 = [];

if exist(filename,'file') ~= 2
    error 'file does not exist';
end

pl2 = internalPL2ReadFileIndex(filename);

end