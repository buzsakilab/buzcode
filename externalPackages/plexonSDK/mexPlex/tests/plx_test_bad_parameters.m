function  [res] = plx_test_bad_parameters()

% trying all functions with invalid file name

res.allTestsPassed = 1;

res=VerifyError('[a,b,c,d,e] = plx_ad(1, 1);', res);
res=VerifyError('[a,b] = plx_ad_chanmap(1);', res);
res=VerifyError('[a,b,c,d] = plx_ad_gap_info(0, 1);', res);
res=VerifyError('[a,b,c] = plx_ad_span(1, 1, -3,-1);', res);
res=VerifyError('[a,b,c] = plx_ad_span_v(1, 1, -3,-1);', res);
res=VerifyError('[a,b,c,d,e] = plx_ad_v(0, 1);', res);
res=VerifyError('[a,b] = plx_adchan_freqs(1);', res);
res=VerifyError('[a,b] = plx_adchan_gains(2);', res);
res=VerifyError('[a,b] = plx_adchan_names(0);', res);
res=VerifyError('[a,b] = plx_adchan_samplecounts(0);', res);

res=VerifyError('[a,b] = plx_chan_filters(1);', res);
res=VerifyError('[a,b] = plx_chan_gains(1);', res);
res=VerifyError('[a,b] = plx_chan_names(2);', res);
res=VerifyError('[a,b] = plx_chan_thresholds(0);', res);
res=VerifyError('[a,b] = plx_chanmap(1);', res);

res=VerifyError('[a,b] = plx_event_chanmap(1);', res);
res=VerifyError('[a,b] = plx_event_names(1);', res);
res=VerifyError('[a,b,c] = plx_event_ts(0, 1);', res);

res=VerifyError('[a,b,c,d]=plx_info(1,0);', res);
res=VerifyError('[a,b,c,d,e,f,g,h,i,j,k,l,m] = plx_information(22);', res);

res=VerifyError('[a,b] = plx_ts(-2, 0, -1);', res);

res=VerifyError('[a,b,c,d] = plx_vt_interpret(1);', res); % here we pass only 1 parameter

res=VerifyError('[a,b,c,d] = plx_waves(0, 1, 1);', res);
res=VerifyError('[a,b,c,d] = plx_waves_v(0, 1, 1);', res);

res=VerifyError('[a,b,c,d] = ddt(1);', res);
res=VerifyError('[a] = ddt_write_v(1, 2, 1, 1, 0);', res);

return;
