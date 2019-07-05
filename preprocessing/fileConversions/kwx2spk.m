function [] =  kwx2spk(basename)



wav = h5read(kwx,'/channel_groups/1/waveforms_filtered');

spktimes = hdf5read(kwik, '/channel_groups/1/spikes/time_samples');



end
