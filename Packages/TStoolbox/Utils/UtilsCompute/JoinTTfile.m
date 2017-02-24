function JoinTTfile(file1, file2, fileout, offset)
  
  
% copy file 1
  
eval(['! cp ' file1 ' ' fileout])




t = LoadTT_tOnly(file2);
n_spikes_tot = length(t);
clear t

spike_offset = 0;
spike_block_size = 50000;


while(spike_offset < n_spikes_tot)
  
  i1 = spike_offset + 1;
  i2 = spike_offset + spike_block_size;
  
  TT_block = LoadTT_chunks(file2, i1, i2);
  
  t2 = Range(TT_block, 'ts') + offset;
  w2 = Data(TT_block);
  
  WriteCheetahTT(fileout, t2, w2, 1);
  spike_offset = spike_offset + spike_block_size;
end
