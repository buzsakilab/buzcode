function [ lfpdown ] = bz_DownsampleLFP( lfp,downsamplefactor )
%[ lfpdown ] = bz_DownsampleLFP( lfp, downsamplefactor )takes a buzcode lfp 
%structure and returns a buzcode lfp structure downsampled by a factor of
%downsamplefactor
%%

lfpdown = lfp;
lfpdown.downsamplefactor = downsamplefactor;
lfpdown.data = downsample(lfpdown.data,lfpdown.downsamplefactor);
lfpdown.samplingRate = lfpdown.samplingRate./lfpdown.downsamplefactor;
lfpdown.timestamps = downsample(lfpdown.timestamps,lfpdown.downsamplefactor);


end

