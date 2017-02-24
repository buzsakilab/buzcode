function  [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(filename)
% plx_information(filename) -- read extended header infromation from a .plx or .pl2 file
%
% [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
% OpenedFileName    - returns the filename (useful if empty string is passed as filename)
% Version -  version code of the plx file format
% Freq -  timestamp frequency for waveform digitization
% Comment - user-entered comment
% Trodalness - 0,1 = single electrode, 2 = stereotrode, 4 = tetrode
% Number of Points Per Wave - number of samples in a spike waveform
% Pre Threshold Points - the sample where the threshold was crossed
% SpikePeakV - peak voltage in mV of the final spike A/D converter
% SpikeADResBits - resolution of the spike A/D converter (usually 12 bits)
% SlowPeakV - peak voltage of mV of the final analog A/D converter
% SlowADResBits - resolution of the analog A/D converter (usually 12 bits)
% Duration - the duration of the file in seconds
% DateTime - date and time string for the file

OpenedFileName = [];
Version = [];
Freq = -1;
Comment = [];
Trodalness = -1;
NPW = -1;
PreTresh = -1;
SpikePeakV = 2500;
SpikeADResBits = 16;
SlowPeakV = 2500;
SlowADResBits = 16;
Duration = -1;
DateTime = [];

if nargin ~= 1
    error 'expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    OpenedFileName = filename;
    Version = pl2.Version.Major*100 + pl2.Version.Minor*10 + pl2.Version.Bugfix;
    Freq = pl2.TimestampFrequency;
    Comment = pl2.CreatorComment;
    if length(Comment) == 0
        Comment = char([]);
    end
    
    DateTime = sprintf('%2d/%2d/%4d %2d:%2d:%2d', pl2.CreatorDateTime.Month+1, pl2.CreatorDateTime.MonthDay, pl2.CreatorDateTime.Year+1900, pl2.CreatorDateTime.Hour, pl2.CreatorDateTime.Minute, pl2.CreatorDateTime.Second);
    Trodalness = pl2.MinimumTrodality;
    if pl2.MaximumTrodality  > pl2.MinimumTrodality 
        Trodalness = [pl2.MinimumTrodality; pl2.MaximumTrodality];
    end
    
    NPW = [];
    PreTresh = [];
    SpikePeakV = [];
    for i=1:numel(pl2.SpikeChannels)
        sps = pl2.SpikeChannels{i}.SamplesPerSpike;
        if length(find(NPW == sps)) == 0
            NPW = [NPW; sps];
        end
        spre = pl2.SpikeChannels{i}.PreThresholdSamples;
        if length(find(PreTresh == spre)) == 0
            PreTresh = [PreTresh; spre];
        end
        vmax = pl2.SpikeChannels{i}.InputVoltageMaximum*1000;
        if length(find(SpikePeakV == vmax)) == 0
            SpikePeakV = [SpikePeakV; vmax];
        end
    end

    SlowPeakV = [];
    for i=1:numel(pl2.AnalogChannels)
        vmax = pl2.AnalogChannels{i}.InputVoltageMaximum*1000;
        if length(find(SlowPeakV == vmax)) == 0
            SlowPeakV = [SlowPeakV; vmax];
        end
    end

    Duration = pl2.DurationOfRecordingTicks/Freq;
    return
end

[OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = mexPlex(13, filename);

end