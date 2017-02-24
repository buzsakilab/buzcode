function [nWritten] = write_plx(filename, ch, freq, npw, n, ts, wave, units)
% write_plx: Write waveform data to a .plx file
%
% [nWritten] = write_plx(filename, ch, freq, npw, n, ts, wave, units)
%
% INPUT:
%   filename - if empty string, will use File Save dialog
%   ch - 1-based channel number that will be written to the .plx file
%   freq - the timestamp frequency in Hz
%   npw - number of points in each waveform
%   n - number of waveforms
%   ts - array of n timestamps (in seconds) 
%   wave - array of waveforms [npw, n], raw a/d values
%   units - array of n unit assignments for each waveform, 0 = unsorted
% OUTPUT:
%   nWritten - number of waveforms successfully written to the output file
%
% This function will write a .plx file that consists of a single channel
% of spike data only. It is not intended to solve the general problem of
% writing out .plx files from Matlab. This code is intended to provide a 
% starting point that can be customized as needed.
%
% Some of the limitations of this function include:
% - There is no facility to write multiple channels, slow a/d data, or event data.
% - Lots of assumptions are made (no comment, non-trodal, prethresh = npw/4,
% the header date is today, etc)
% - The count arrays in the header (TSCounts, WFCounts, EVCounts) are not
% filled (this is fixable with PlexUtil)
% - the gain for the channel is set to 1, which will result in the 
% reconstructed voltages for the spikes to be different.
% - unit templates are not written out, they can be reconstructed with OFS
% - long timestamps are not supported
% 
% Note that the PlexUtil program from Plexon can merge the output files 
% from this routine together to form multi-channel files.
%
%


if(length(filename) == 0)
   [fname, pathname] = uigetfile('*.plx', 'Select a plx file');
	filename = strcat(pathname, fname);
end

pad256(1:256) = uint8(0);

% create the file and write the file header

plx_id = fopen(filename, 'w');
fwrite(plx_id, 1480936528, 'integer*4');    % 'PLEX' magic code
fwrite(plx_id, 101, 'integer*4');           % the version no.
fwrite(plx_id, pad256(1:128), 'char');      % placeholder for comment
fwrite(plx_id, freq, 'integer*4');          % timestamp frequency
fwrite(plx_id, 1, 'integer*4');             % no. of DSP channels
fwrite(plx_id, 0, 'integer*4');             % no. of event channels
fwrite(plx_id, 0, 'integer*4');             % no. of A/D (slow-wave) channels
fwrite(plx_id, npw, 'integer*4');           % no. points per waveform
fwrite(plx_id, npw/4, 'integer*4');         % (fake) no. pre-threshold points
[YR, MO, DA, HR, MI, SC] = datevec(now);    % current date & time
fwrite(plx_id, YR, 'integer*4');            % year
fwrite(plx_id, MO, 'integer*4');            % month
fwrite(plx_id, DA, 'integer*4');            % day
fwrite(plx_id, HR, 'integer*4');            % hour
fwrite(plx_id, MI, 'integer*4');            % minute
fwrite(plx_id, SC, 'integer*4');            % second
fwrite(plx_id, 0, 'integer*4');             % fast read (reserved)
fwrite(plx_id, freq, 'integer*4');          % waveform frequency
fwrite(plx_id, max(ts)*freq, 'double');     % last timestamp
fwrite(plx_id, pad256(1:56), 'char');       % should make 256 bytes

% now the count arrays (with counts of zero)

for i = 1:40
  fwrite(plx_id, pad256(1:130), 'char');    % first 20 are TSCounts, next 20 are WFCounts
end
for i = 1:8
  fwrite(plx_id, pad256(1:256), 'char');    % all of these make up EVCounts
end

% now the single PL_ChanHeader

% assume simple channel names
[DSPname, errmsg] = sprintf('DSP%d',ch);
[SIGname, errmsg] = sprintf('SIG%d',ch);

fwrite(plx_id, DSPname, 'char');
fwrite(plx_id, pad256(1:32-length(DSPname)));
fwrite(plx_id, SIGname, 'char');
fwrite(plx_id, pad256(1:32-length(SIGname)));
fwrite(plx_id, ch, 'integer*4');            % DSP channel number
fwrite(plx_id, 0, 'integer*4');             % waveform rate limit (not used)
fwrite(plx_id, ch, 'integer*4');            % SIG associated channel number
fwrite(plx_id, ch, 'integer*4');            % SIG reference  channel number
fwrite(plx_id, 1, 'integer*4');             % dummy for gain
fwrite(plx_id, 0, 'integer*4');             % filter off
fwrite(plx_id, -12, 'integer*4');           % (fake) detection threshold value
fwrite(plx_id, 0, 'integer*4');             % sorting method (dummy)
fwrite(plx_id, 0, 'integer*4');             % number of sorted units
for i = 1:10
  fwrite(plx_id, pad256(1:64), 'char');     % filler for templates (5 * 64 * short)
end
fwrite(plx_id, pad256(1:20), 'char');       % template fit (5 * int)
fwrite(plx_id, npw, 'integer*4');           % sort width (template only)
fwrite(plx_id, pad256(1:80), 'char');       % boxes (5 * 2 * 4 * short)
fwrite(plx_id, 0, 'integer*4');             % beginning of sorting window
fwrite(plx_id, pad256(1:128), 'char');      % comment
fwrite(plx_id, pad256(1:44), 'char');       % padding

% now the spike waveforms, each preceded by a PL_DataBlockHeader

for ispike = 1:n

  fwrite(plx_id, 1, 'integer*2');           % type: 1 = spike
  fwrite(plx_id, 0, 'integer*2');           % upper byte of 5-byte timestamp
  fwrite(plx_id, ts(ispike)*freq, 'integer*4');  % lower 4 bytes
  fwrite(plx_id, ch, 'integer*2');          % channel number
  fwrite(plx_id, units(ispike), 'integer*2');  % unit no. (0 = unsorted)
  fwrite(plx_id, 1, 'integer*2');           % no. of waveforms = 1
  fwrite(plx_id, npw, 'integer*2');         % no. of samples per waveform

  fwrite(plx_id, wave(ispike, 1:npw), 'integer*2');
   
end

nWritten = n;

fclose(plx_id);
