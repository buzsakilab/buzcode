function [] = plx_spike_info(filename)
% plx_spike_info(filename): prints .plx file spike channel info
%
% plx_spike_info(filename)
%
% INPUT:
%   filename - .plx file name. Will use file open dialog if filename is empty string
%
% OUTPUT:
%   (none)

if nargin ~= 1
    error 'expected 1 input argument';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    PL2Print(pl2.SpikeChannels);
    return
end

[n, names] = plx_chan_names(filename);
[n, raw] = plx_chanmap(filename);
[tscounts, wfcounts, evcounts, contcounts] = plx_info(filename, 1);
fprintf(' %18.18s', 'Channel Number');
fprintf(' %18.18s', 'Name');
fprintf(' %18.18s', 'Total Waveforms');
fprintf(' %18.18s', 'Unsorted');
fprintf(' %18.18s', 'Unit a');
fprintf(' %18.18s', 'Unit b');
fprintf(' %18.18s', 'Unit c');
fprintf(' %18.18s', 'Unit d');
fprintf('\n');
    
for i=1:n
    fprintf(' %18d', raw(i));
    fprintf(' %18.18s', names(i,:));
    fprintf(' %18d', sum(wfcounts(:,i+1)));
    fprintf(' %18d', wfcounts(1,i+1));
    fprintf(' %18d', wfcounts(2,i+1));
    fprintf(' %18d', wfcounts(3,i+1));
    fprintf(' %18d', wfcounts(4,i+1));
    fprintf(' %18d', wfcounts(5,i+1));
    fprintf(' %18d', wfcounts(6,i+1));
    fprintf(' %18d', wfcounts(7,i+1));
    fprintf(' %18d', wfcounts(8,i+1));
    fprintf(' %18d', wfcounts(9,i+1));
    fprintf(' %18d', wfcounts(10,i+1));
    fprintf(' %18d', wfcounts(11,i+1));
    fprintf('\n');
end

end