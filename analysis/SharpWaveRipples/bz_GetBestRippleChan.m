function  [chan] = bz_GetBestRippleChan(lfp)
%[chan] = bz_GetBestRippleChan(lfp)
%eventually this will detect which lfp channel has the highest SNR for the
% ripple componenent of SPWR events....

[b a]=butter(4,[140/625 180/625],'bandpass');

for i=1:length(lfp.channels)
    filt = FiltFiltM(b,a,single(lfp.data(:,i)));
    pow = fastrms(filt,15);    
    mRipple(i) = mean(pow);
    meRipple(i) = median(pow);
    mmRippleRatio(i) = mRipple(i)./meRipple(i);
end

mmRippleRatio(mRipple<1) = 0;
mmRippleRatio(meRipple<1) = 0;

[minVal loc] = max(mmRippleRatio);
chan = lfp.channels(loc);
end