% this script tests the import routines
% IT IS NOT USEFUL IN GENERAL, as it makes assumptions about what data is
% in the plx file.
disp('This sample .m script requires customization to work properly!');

% Open a plx file
% this will bring up the file-open dialog
StartingFileName = '';
%StartingFileName = 'C:\PlexonData\CM_Quickstart.plx';
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

% get some counts
[tscounts, wfcounts, evcounts, slowcounts] = plx_info(OpenedFileName,1);

% tscounts, wfcounts are indexed by (channel+1,unit+1)
% tscounts(:,ch+1) is the per-unit counts for channel ch
% sum( tscounts(:,ch+1) ) is the total wfs for channel ch (all units)
% [nunits, nchannels+1] = size( tscounts )
% To get number of nonzero units/channels, use nnz() function

% get some strobed timestamps
[nev, evts, evsv] = plx_event_ts(OpenedFileName, 257);
if nev > 0 
    [nCoords, nDim, nVTMode, c] = plx_vt_interpret(evts, evsv);
    disp('VT data is interpreted.');
end

% get some timestamps for channel 1 unit a
[nts, ts] = plx_ts(OpenedFileName, 1, 1);
[nwf, npw, tswf, waves] = plx_waves(OpenedFileName, 1, 1);
if nwf > 0
    plot( waves(1,1:NPW));
end

% get some other info about the spike channels
[nspkfilters,spk_filters] = plx_chan_filters(OpenedFileName);
[nspkgains,spk_gains] = plx_chan_gains(OpenedFileName);
[nspkthresh,spk_threshs] = plx_chan_thresholds(OpenedFileName);

% get some strobed timestamps
%[nev, evts, evsv] = plx_event_ts(OpenedFileName, 257);
% get some non-strobed timestamps
[nev, evts, evsv] = plx_event_ts(OpenedFileName, 2);

% get some a/d data
[adfreq, nad, tsad, fnad, ad] = plx_ad(OpenedFileName, 1);
if nad > 0
    plot( ad );
end

% get just a span of a/d data
[adfreq, nadspan, adspan] = plx_ad_span(OpenedFileName, 1, 10,100);

[nadfreqs,adfreqs] = plx_adchan_freqs(OpenedFileName);
[nadgains,adgains] = plx_adchan_gains(OpenedFileName);

