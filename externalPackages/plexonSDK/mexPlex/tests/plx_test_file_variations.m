function  [res] = plx_test_file_variations()

% trying functions with files with some header and strobed parameters

res.allTestsPassed = 0;

res=VerifyError('[a,b,c,d] = plx_info(''ts_freq_zero.plx\'',1);', res);
[a,b,c,d] = plx_info('strobed_negative.plx',1);
[a,b,c,d] = plx_info('waveform_freq_zero.plx',1);

res.allTestsPassed = 1;

return;
