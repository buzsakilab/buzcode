% this script tests the import routines
% IT IS NOT USEFUL IN GENERAL, as it makes assumptions about what data is
% in the plx file.
disp('This sample .m script requires customization to work properly!');

% Open a plx file
% this will bring up the file-open dialog
%StartingFileName = '';
StartingFileName = 'C:\PlexonData\CM_Quickstart.plx';
%StartingFileName = 'c:\plexondata\NSSample.plx';
[OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreThresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(StartingFileName);

disp(['Opened File Name: ' OpenedFileName]);
disp(['Version: ' num2str(Version)]);
disp(['Frequency : ' num2str(Freq)]);
disp(['Comment : ' Comment]);
disp(['Date/Time : ' DateTime]);
disp(['Duration : ' num2str(Duration)]);
disp(['Num Pts Per Wave : ' num2str(NPW)]);
disp(['Num Pts Pre-Threshold : ' num2str(PreThresh)]);
% some of the information is only filled if the plx file version is >102
if ( Version > 102 )
    if ( Trodalness < 2 )
        disp('Data type : Single Electrode');
    elseif ( Trodalness == 2 )
        disp('Data type : Stereotrode');
    elseif ( Trodalness == 4 )
        disp('Data type : Tetrode');
    else
        disp('Data type : Unknown');
    end
        
    disp(['Spike Peak Voltage (mV) : ' num2str(SpikePeakV)]);
    disp(['Spike A/D Resolution (bits) : ' num2str(SpikeADResBits)]);
    disp(['Slow A/D Peak Voltage (mV) : ' num2str(SlowPeakV)]);
    disp(['Slow A/D Resolution (bits) : ' num2str(SlowADResBits)]);
end   


% get the unsorted waveforms
channel = 1;
[n, npw, ts, wave] = plx_waves(OpenedFileName, channel, 0);

% set up a fake units array, 'sort' the spikes into unit 2
units(1:size(ts)) = 2;

% make an output file name
OutputFileName = 'C:\PlexonData\plx_write_test\plx_write_out_ch1.plx';

% write it out
[nWritten] = write_plx(OutputFileName, channel, Freq, npw, n, ts, wave, units)


% do another channel
channel = 2;
[n, npw, ts, wave] = plx_waves(OpenedFileName, channel, 0);

% set up a fake units array, 'sort' the spikes into unit 1
units(1:size(ts)) = 1;

% make an output file name
OutputFileName = 'C:\PlexonData\plx_write_test\plx_write_out_ch2.plx';

% write it out
[nWritten] = write_plx(OutputFileName, channel, Freq, npw, n, ts, wave, units)
