function writePLXFile(filename, data, h)
% writePLXFile - Write data to a PLX file.
%
% writePLXFile(filename, data, headers)
%
% writePLXFile is designed to generate a new PLX file using spiking data
% and header information. The header data can be taken from an existing PLX
% file using the companion function readPLXHeaders, or reasonable default
% values will be created by checkPLXData. This function relies upon
% checkPXLData to verify that the data is in the correct format and that
% the headers match with the data before writing.
%
% A typical use case is to read a data file using the libraries provided by
% Plexon, manipulate the data using custom written MATLAB functions, then
% write the manipulated data back to a new PLX file, using the headers read
% from the original PLX file. checkPLXData will adjust the appropriate
% fields in the headers to account for changes in the data (such as the
% number of samples per spike) before writing the data file.
%
% WARNING: Continuous (slow channel) and event data is not yet supported by
% this function. All related data is thrown away and all related fields are
% zeroed out.
%
% Input:
%   filename - Filename of newly created PLX file.
%   data - A structure containing the spiking data
%     data.ts    - [n x 1] array of timestamps (seconds)
%     data.chans - [n x 1] array of channel assignments (1 based)
%     data.units - [n x 1] array of unit assignments (0 = unsorted)
%     data.waves - [n x npw] matrix of waveforms
%   headers - PLX file header information (optional). See readPLXHeaders.
%
% See also readPLXHeaders, checkPLXData.
%
% AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% Copyright (c) 2011, Benjamin Kraus
% $Id$

if(isempty(filename))
   [fname, pathname] = uiputfile('*.plx', 'Select a plx file');
	filename = strcat(pathname, fname);

end

if(nargin == 2)
    h = checkPLXData(data);
elseif(nargin == 3)
    h = checkPLXData(data, h);
else
    error('Two or three arguments expected');
end

nDSP = h.numDSPChannels;
npw = h.numPointsWave;

% create the file and write the file header
fid = fopen(filename, 'w');
if(fid == -1); error('writePLXFile:FileError','Error opening file for writing'); end
c = onCleanup(@()fclose(fid));

fwrite(fid, 'PLEX', 'char');                  % 'PLEX' magic code
fwrite(fid, h.version, 'int32');              % the version no.
fwrite(fid, padstr(h.comment, 128), 'char');  % placeholder for comment
fwrite(fid, h.ADFrequency, 'int32');          % timestamp frequency
fwrite(fid, nDSP, 'int32');                   % no. of DSP channels
fwrite(fid, 0, 'int32');                      % no. of event channels
fwrite(fid, 0, 'int32');                      % no. of A/D (slow-wave) channels
fwrite(fid, npw, 'int32');                    % no. points per waveform
fwrite(fid, h.numPointsPreThr, 'int32');      % (fake) no. pre-threshold points

[YR, MO, DA, HR, MI, SC] = datevec(h.date);
fwrite(fid, YR, 'int32');  % year
fwrite(fid, MO, 'int32'); % month
fwrite(fid, DA, 'int32'); % day
fwrite(fid, HR, 'int32'); % hour
fwrite(fid, MI, 'int32'); % minute
fwrite(fid, SC, 'int32'); % second

fwrite(fid, h.fastread, 'int32');       % fast read (reserved)
fwrite(fid, h.waveformfreq, 'int32');   % waveform frequency
fwrite(fid, h.lasttimestamp, 'double'); % last timestamp

fwrite(fid, h.trodalness, 'uint8');           % trodalness
fwrite(fid, h.datatrodalness, 'uint8');       % data trodalness
fwrite(fid, h.bitsperspikesample, 'uint8');   % bits per spike sample
fwrite(fid, h.bitsperslowsample, 'uint8');    % bits per slow sample
fwrite(fid, h.spikeMaxMagnitudeMV, 'ushort'); % spike max magnitude (mV)
fwrite(fid, h.slowMaxMagnitudeMV, 'ushort');  % slow max magnitude (mV)

fwrite(fid, h.spikePreAmpGain, 'ushort');     % spike pre-amp gain

fwrite(fid, padstr(h.acquiringsoftware,18), 'char');  % acquiring software
fwrite(fid, padstr(h.processingsoftware,18), 'char'); % processing software

fwrite(fid, padstr('',10),'char'); % should total 256 bytes

% now the count arrays (with counts of zero)
fwrite(fid, h.tscounts(1:5,1:130),'int32');
fwrite(fid, h.wfcounts(1:5,1:130),'int32');
fwrite(fid, h.evcounts(1:512),'int32');


for ch = 1:nDSP;
    fwrite(fid, padstr(h.chans(ch).name,32), 'char');
    fwrite(fid, padstr(h.chans(ch).signame,32), 'char');
    fwrite(fid, h.chans(ch).channel, 'int32');              % DSP channel number
    fwrite(fid, h.chans(ch).wfrate, 'int32');               % waveform rate limit (not used)
    fwrite(fid, h.chans(ch).sig, 'int32');                  % SIG associated channel number
    fwrite(fid, h.chans(ch).ref, 'int32');                  % SIG reference  channel number
    fwrite(fid, h.chans(ch).gain, 'int32');                 % dummy for gain
    fwrite(fid, h.chans(ch).filter, 'int32');               % filter off
    fwrite(fid, h.chans(ch).threshold, 'int32');            % (fake) detection threshold value
    fwrite(fid, h.chans(ch).method, 'int32');               % sorting method (dummy)
    fwrite(fid, h.chans(ch).nunits, 'int32');               % number of sorted units
    fwrite(fid, h.chans(ch).template, 'short');
    fwrite(fid, h.chans(ch).fit, 'int32');                   % template fit (5 * int)
    fwrite(fid, h.chans(ch).sortwidth, 'int32');            % sort width (template only)
    fwrite(fid, h.chans(ch).boxes, 'short');                 % boxes (5 * 2 * 4 * short)
    fwrite(fid, h.chans(ch).sortbeg, 'int32');              % beginning of sorting window
    fwrite(fid, padstr(h.chans(ch).comment,128), 'char');   % comment
    fwrite(fid, padstr(h.chans(ch).srcid,1), 'uchar');
    fwrite(fid, char(0), 'uchar');
    fwrite(fid, h.chans(ch).chanid,'ushort');
    fwrite(fid, zeros(10,1), 'int32');
end

[~,ind] = sortrows([data.ts, data.chans],[1 2]);
data.ts = data.ts(ind,:);
data.chans = data.chans(ind,:);
data.units = data.units(ind,:);
data.waves = data.waves(ind,:);

ts = round(data.ts*h.ADFrequency);
h = waitbar(0,sprintf('Writing Spiking Data (%02.0f%%)',0));
nspks = size(data.ts,1);
updateintv = floor(nspks/100);
for ii = 1:nspks;
    fwrite(fid, 1, 'short');           % type: 1 = spike
    fwrite(fid, 0, 'short');           % upper byte of 5-byte timestamp
    fwrite(fid, ts(ii), 'int32');  % lower 4 bytes
    fwrite(fid, data.chans(ii), 'short');          % channel number
    fwrite(fid, data.units(ii), 'short');  % unit no. (0 = unsorted)
    fwrite(fid, 1, 'short');           % no. of waveforms = 1
    fwrite(fid, npw, 'short');         % no. of samples per waveform
    fwrite(fid, data.waves(ii, 1:npw), 'short');
    if(mod(ii,updateintv)==0);
        x = ii/nspks;
        waitbar(x,h,sprintf('Writing Spiking Data (%02.0f%%)',x*100));
    end
end
close(h);

end

function strout = padstr(strin, len)
    strin = strin(:)';
    if(length(strin)>len); strout = strin(1:len);
    else strout = [strin char(zeros(1, len-length(strin)))];
    end
end