function  [plx] = get_all_from_plx(fileName)

plx = [];
if (nargin == 0 | isempty(fileName))
   [fname, pathname] = uigetfile('*.plx', 'Select a .plx file');
   if isequal(fname,0)
     error 'No file was selected'
     return
   end
   fileName = fullfile(pathname, fname);
end

[plx.tscounts, plx.wfcounts, plx.evcounts, plx.slowcounts] = plx_info(fileName,1);
[OpenedFilePath, plx.Version, plx.Freq, plx.Comment, plx.Trodalness, plx.NPW, plx.PreThresh, plx.SpikePeakV, plx.SpikeADResBits, plx.SlowPeakV, plx.SlowADResBits, plx.Duration, plx.DateTime] = plx_information(fileName);
[dirstr, name, ext] = fileparts(OpenedFilePath);
plx.FileName = [name ext];

% try valid and invalid channels and units
% max unit is 26
% max spike channel number is 128
for iunit = -1:30  
    for ich = -1:130
        [plx.ts.n{iunit+2,ich+2}, plx.ts.ts{iunit+2,ich+2}] = plx_ts(fileName, ich , iunit );
        [plx.wf.n{iunit+2,ich+2}, plx.wf.npw{iunit+2,ich+2}, plx.wf.ts{iunit+2,ich+2}, plx.wf.wf{iunit+2,ich+2}] = plx_waves(fileName, ich , iunit );
        [plx.wf_v.n{iunit+2,ich+2}, plx.wf_v.npw{iunit+2,ich+2}, plx.wf_v.ts{iunit+2,ich+2}, plx.wf_v.wf{iunit+2,ich+2}] = plx_waves_v(fileName, ich , iunit );
     end
end

[nspk,plx.spk_filters] = plx_chan_filters(fileName);
[nspk,plx.spk_gains] = plx_chan_gains(fileName);
[nspk,plx.spk_threshs] = plx_chan_thresholds(fileName);
[nspk,plx.spk_names] = plx_chan_names(fileName);

% try valid and invalid channels
% ona 128-channel system, there could be 128*3 a/d channels
for ich = -1:400
	[plx.ad_gi.freq{ich+2}, plx.ad_gi.nad{ich+2}, plx.ad_gi.ts{ich+2}, plx.ad_gi.frag{ich+2}] = plx_ad_gap_info(fileName, ich);
	[plx.ad.freq{ich+2}, plx.ad.nad{ich+2}, plx.ad.ts{ich+2}, plx.ad.frag{ich+2}, plx.ad.val{ich+2}] = plx_ad(fileName, ich);
	[plx.ad_v.freq{ich+2}, plx.ad_v.nad{ich+2}, plx.ad_v.ts{ich+2}, plx.ad_v.frag{ich+2}, plx.ad_v.val{ich+2}] = plx_ad_v(fileName, ich);
	[plx.ad_span.freq{ich+2}, plx.ad_span.nad{ich+2}, plx.ad_span.val{ich+2}] = plx_ad_span(fileName, ich, 10,20);
	[plx.ad_span_v.freq{ich+2}, plx.ad_span_v.nad{ich+2} plx.ad_span_v.val{ich+2}] = plx_ad_span_v(fileName, ich, 10,20);
    % try invalid start/end combinations
	[plx.ad_span_inv1.freq{ich+2}, plx.ad_span_inv1.nad{ich+2} plx.ad_span_inv1.val{ich+2}] = plx_ad_span_v(fileName, ich, 10,-5);
    [plx.ad_span_inv2.freq{ich+2}, plx.ad_span_inv2.nad{ich+2} plx.ad_span_inv2.val{ich+2}] = plx_ad_span_v(fileName, ich, 10,9);
    [plx.ad_span_inv3.freq{ich+2}, plx.ad_span_inv3.nad{ich+2} plx.ad_span_inv3.val{ich+2}] = plx_ad_span_v(fileName, ich, -2,3);
end		

[nad,plx.adfreqs] = plx_adchan_freqs(fileName);
[nad,plx.adgains] = plx_adchan_gains(fileName);
[nad,plx.adnames] = plx_adchan_names(fileName);   

% the events
[u,plx.evnames] = plx_event_names(fileName);
[u,plx.evchans] = plx_event_chanmap(fileName);
for iev = -1:310
	[plx.nevs{iev+2}, plx.tsevs{iev+2}, plx.svdummy{iev+2}] = plx_event_ts(fileName, iev);
	if ( iev > 0 && iev <= length(plx.evchans) && plx.evcounts(iev) > 0 )
	    evch = plx.evchans(iev);
	    if ( evch == 257 )
			[plx.strobed.n, plx.strobed.ts, plx.strobed.v] = plx_event_ts(fileName, evch); 
            [plx.nCoords, plx.nDim, plx.nVTMode, plx.c] = plx_vt_interpret(plx.strobed.ts, plx.strobed.v);
        end
    end
end

plx_close(fileName);

% try to get a/d channel after closing
for ich = -1:400
	[plx.ad2.freq{ich+2}, plx.ad2.nad{ich+2}, plx.ad2.ts{ich+2}, plx.ad2.frag{ich+2}, plx.ad2.val{ich+2}] = plx_ad(fileName, ich);
end
plx_close('');