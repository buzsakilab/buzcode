function  [tscounts, wfcounts, evcounts, contcounts] = plx_info(filename, fullread)
% plx_info(filename, fullread) -- read and display .plx or .pl2 file info
%
% [tscounts, wfcounts, evcounts, contcounts] = plx_info(filename, fullread)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   fullread - if 0, reads only the file header
%              if 1, reads the entire file
%               for .pl2 files, this parameter is ignored
%
% OUTPUT:
%   tscounts - 2-dimensional array of timestamp counts for each unit
%      tscounts(i, j) is the number of timestamps for channel j-1, unit i
%                                (see comment below)
%   wfcounts - 2-dimensional array of waveform counts for each unit
%     wfcounts(i, j) is the number of waveforms for channel j-1, unit i
%                                (see comment below)
%   evcounts - 1-dimensional array of external event counts
%     evcounts(i) is the number of events for event channel i
%
%   contcounts - 1-dimensional array of sample counts for continuous channels
%     contcounts(i) is the number of continuous for slow channel i-1
%
% Note that for tscounts, wfcounts, the unit,channel indices i,j are off by one. 
% That is, for channels, the count for channel n is at index n+1, and for units,
%  index 1 is unsorted, 2 = unit a, 3 = unit b, etc
% The dimensions of the tscounts and wfcounts arrays are
%   (NChan+1) x (MaxUnits+1)
% where NChan is the number of spike channel headers in the plx file, and
% MaxUnits is 4 if fullread is 0, or 26 if fullread is 1. This is because
% the header of a .plx file can only accomodate 4 units, but doing a
% fullread on the file may show that there are actually up to 26 units
% present in the file. Likewise, NChan will have a maximum of 128 channels
% if fullread is 0.
% The dimension of the evcounts and contcounts arrays is the number of event
% and continuous (slow) channels. 
% The counts for slow channel 0 is at contcounts(1)

tscounts = [];
wfcounts = [];
evcounts = [];
contcounts = [];

if nargin ~= 2
    error 'expected 2 input arguments';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    pl2 = PL2GetFileIndex(filename);
    numSpikeChannels = numel(pl2.SpikeChannels);
    % pl2 files support up to 256 units, but we limit to 
    % 26 sorted plus 1 unsorted to be compatible with plx_info
    tscounts = zeros(27,numSpikeChannels+1);
    for i=1:numSpikeChannels
        tscounts(:,i+1) = pl2.SpikeChannels{i}.UnitCounts(1:27);
    end
    wfcounts = tscounts;

    numAnalogChannels = numel(pl2.AnalogChannels);
    contcounts = zeros(1,numAnalogChannels);
    for i=1:numAnalogChannels
        contcounts(1,i) = pl2.AnalogChannels{i}.NumValues;
    end

    evcounts = pl2.EventCounts;
    return
end

[tscounts, wfcounts, evcounts, contcounts] = mexPlex(4, filename, fullread);

end