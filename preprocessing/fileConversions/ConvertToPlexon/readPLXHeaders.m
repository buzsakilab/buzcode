function [headers, fid] = readPLXHeaders(filename, fullread)
% readPLXHeaders - Read file headers and channel information from PLX file.
%
% [headers, fid] = readPLXHeaders(filename)
%
% readPLXHeaders will read the header data from a PLX file and return a
% structure with fields corresponding to each piece of header data.
%
% Input:
%   filename - Filename of a PLX file. A valid file-ID can be provided
%              instead of a filename. In this case, the file-ID will not be
%              closed after reading.
%   fullread - Verify the contents of the file header by reading the entire
%              file. This effects the values of headers.tscounts,
%              headers.wfcounts, headers.evcounts, and headers.slowcounts.
%              It will also verify that headers.chans(n).nunits is
%              accurate, and that headers.numPointsWave is accurate.
%              (default = false)
%
% Output:
%   headers - PLX file header information.
%   fid     - File-ID of the 
%
% WARNING: Plexon didn't have much foresight for the number of
% units/channel, and only left enough space in the headers for 4
% units/channel. Therefore, the data stored in the headers (and output as
% headers.tscounts and headers.wfcounts) can't account for more than 4
% units/channel. The only way to get an accurate count of the number of
% timestamps, waveforms, and events is to read the entire data file. If you
% call this function with the 'fullread' argument set false (default) it
% will only read the values provided in the header. If you set 'fullread'
% to true, it will read the entire file to verify the information in the
% header. This is similar to the behavior of the Plexon function called
% plx_info.
%
% Note: This function behaves slightly different than plx_info. In
% particular:
%     tscounts   = 27 x m array. The entry at tscounts(c,u) refers to
%                  channel c and unit u-1. Channel numbers start at 1.
%                  Unit 0 (in row 1) is unsorted waves.
%     wfcounts   = 27 x m array. The entry at wfcounts(c,u) refers to
%                  channel c and unit u-1. Channel numbers start at 1.
%                  Unit 0 (in row 1) is unsorted waves.
%     evcounts   = 1 x m array. The entry at evcounts(c) refers to
%                  channel c. Channel numbers start at 1. m is the event
%                  with the largest channel number (usually 259).
%     slowcounts = 1 x m+1 array. The entry at slowcounts(c) refers to
%                  channel c+1. Channel numbers start at 0. m is one
%                  greater than the largest slow channel number.
%
% See also writePLXFile, checkPLXData, plx_info, plx_information.
%
% AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% Copyright (c) 2011, Benjamin Kraus
% $Id$

DEBUG = true;

if (nargin == 0 || isempty(filename))
    FilterSpec = {'*.plx', 'Plexon PLX File (*.plx)';
                  '*', 'All Files'};
    [fname, pathname] = uigetfile(FilterSpec, 'Select a Plexon PLX file');        
    if(fname == 0); headers = struct([]); return; end
    filename = strcat(pathname, fname);
end

if(~ischar(filename));
    % If the filename input is not a character string, check to see if it
    % is a valid file ID instead. We can do this by calling 'frewind'. If
    % 'frewind' succeeds, then we know we have a valid file ID.
    fid = filename;
    frewind(fid);
else
    fid = fopen(filename, 'r');
    if(fid ~= -1 && nargout < 2); c = onCleanup(@()fclose(fid)); end
end

if(fid == -1); error('readPLXHeaders:FileError','Error opening file'); end

if(nargin < 2); fullread = false; shortread = false;
elseif(islogical(fullread)); shortread = false;
elseif(ischar(fullread))
    switch lower(fullread)
        case 'true';   fullread = true;  shortread = false;
        case 'false';  fullread = false; shortread = false;
        case 'full';   fullread = true;  shortread = false;
        case 'short';  fullread = false; shortread = true;
        otherwise;     fullread = false; shortread = false;
    end
end

% Determine file size:
fseek(fid,0,'eof');
fsize = ftell(fid);
frewind(fid);

% File header is 256 bytes

magic = fread(fid, 4, '*char')';
if(~strcmp(magic,'PLEX'))
    error('readPLXHeaders:InvalidPLXFile','This file is not a valid .plx file');
end

headers.version = fread(fid, 1, 'int32');
headers.comment = deblank(fread(fid, 128, '*char')');
headers.ADFrequency = fread(fid, 1, 'int32');
headers.numDSPChannels = fread(fid, 1, 'int32');
headers.numEventChannels = fread(fid, 1, 'int32');
headers.numSlowChannels = fread(fid, 1, 'int32');
headers.numPointsWave = fread(fid, 1, 'int32');
headers.numPointsPreThr = fread(fid, 1, 'int32');

YR = fread(fid, 1, 'int32');
MO = fread(fid, 1, 'int32');
DA = fread(fid, 1, 'int32');
HR = fread(fid, 1, 'int32');
MI = fread(fid, 1, 'int32');
SC = fread(fid, 1, 'int32');
headers.date = datenum([YR, MO, DA, HR, MI, SC]);

headers.fastread = fread(fid, 1, 'int32');
headers.waveformfreq = fread(fid, 1, 'int32');
headers.lasttimestamp = fread(fid, 1, 'double');

if(headers.version >= 103)
    headers.trodalness = fread(fid, 1, 'uint8');
    headers.datatrodalness = fread(fid, 1, 'uint8');
    headers.bitsperspikesample = fread(fid, 1, 'uint8');
    headers.bitsperslowsample = fread(fid, 1, 'uint8');
    headers.spikeMaxMagnitudeMV = fread(fid, 1, 'ushort');
    headers.slowMaxMagnitudeMV = fread(fid, 1, 'ushort');
else
    headers.trodalness = 1;
    headers.datatrodalness = 1;
    headers.bitsperspikesample = 12;
    headers.bitsperslowsample = 12;
    headers.spikeMaxMagnitudeMV = 3000;
    headers.slowMaxMagnitudeMV = 5000;
    fread(fid, 8, 'uint8');
end

if(headers.version >= 105)
    headers.spikePreAmpGain = fread(fid, 1, 'ushort');
else
    headers.spikePreAmyGain = 1000;
    fread(fid, 1, 'ushort');
end

if(headers.version >= 106)
    headers.acquiringsoftware = deblank(fread(fid, 18, '*char')');
    headers.processingsoftware = deblank(fread(fid, 18, '*char')');
    fread(fid, 10, 'uint8');
else
    headers.acquiringsoftware = '';
    headers.processingsoftware = '';
    fread(fid, 46, 'char');
end

% Current position = 256

% Timestamps, waveforms, event counts total 7248 bytes.

tscounts = fread(fid, [5, 130], 'int32');
wfcounts = fread(fid, [5, 130], 'int32');
evcounts = fread(fid, [1, 512], 'int32');

headers.tscounts = tscounts(:,2:min(headers.numDSPChannels+1,130));
headers.wfcounts = wfcounts(:,2:min(headers.numDSPChannels+1,130));
headers.evcounts = evcounts(2:300);
headers.slowcounts = evcounts(301:min(300+headers.numSlowChannels,512));

% Current position = 7504

% If instructed to do a 'short' read, stop here.
if(shortread); return; end

% Each channel header is 1020 bytes.
headers.chans = struct([]);
for ii = 1:headers.numDSPChannels;
    headers.chans(ii,1).name = deblank(fread(fid, 32, '*char')');
    headers.chans(ii,1).signame = deblank(fread(fid, 32, '*char')');
    headers.chans(ii,1).channel = fread(fid, 1, 'int32');
    headers.chans(ii,1).wfrate = fread(fid, 1, 'int32');
    headers.chans(ii,1).sig = fread(fid, 1, 'int32');
    headers.chans(ii,1).ref = fread(fid, 1, 'int32');
    headers.chans(ii,1).gain = fread(fid, 1, 'int32');
    headers.chans(ii,1).filter = fread(fid, 1, 'int32');
    headers.chans(ii,1).threshold = fread(fid, 1, 'int32');
    headers.chans(ii,1).method = fread(fid, 1, 'int32');
    headers.chans(ii,1).nunits = fread(fid, 1, 'int32');
    headers.chans(ii,1).template = fread(fid, [5, 64], 'short');
    headers.chans(ii,1).fit = fread(fid, 5, 'int32');
    headers.chans(ii,1).sortwidth = fread(fid, 1, 'int32');
    headers.chans(ii,1).boxes = zeros(5, 2, 4);
    headers.chans(ii,1).boxes(:) = fread(fid, 5*2*4, 'short');
    headers.chans(ii,1).sortbeg = fread(fid, 1, 'int32');
    if(headers.version>=105)
        headers.chans(ii,1).comment = deblank(fread(fid, 128, '*char')');
    else
        headers.chans(ii,1).comment = '';
        fread(fid, 128, 'char');
    end
    
    if(headers.version>=106)
        headers.chans(ii,1).srcid = fread(fid, 1, 'uchar');
        fread(fid, 1, 'uchar');
        headers.chans(ii,1).chanid = fread(fid, 1, 'ushort');
        fread(fid, 10, 'int32');
    else fread(fid, 11, 'int32');
    end
end

% Current position = 7504 + 1020*nchans

% Each event channel header is 296 bytes
headers.evchans = struct([]);
for ii = 1:headers.numEventChannels;
    headers.evchans(ii,1).name = deblank(fread(fid, 32, '*char')');
    headers.evchans(ii,1).num = fread(fid, 1, 'int32');

    if(headers.version>=105)
        headers.evchans(ii,1).comment = deblank(fread(fid, 128, '*char')');
    else
        headers.evchans(ii,1).comment = '';
        fread(fid, 128, 'char');
    end
    
    if(headers.version>=106)
        headers.evchans(ii,1).srcid = fread(fid, 1, 'uchar');
        fread(fid, 1, 'uchar');
        headers.evchans(ii,1).chanid = fread(fid, 1, 'ushort');
        fread(fid, 32, 'int32');
    else fread(fid, 33, 'int32');
    end
end

% Current position = 7504 + 1020*nchans + 296*nevchans

% Each slow channel header is 296 bytes
headers.slowchans = struct([]);
for ii = 1:headers.numSlowChannels;
    headers.slowchans(ii,1).name = deblank(fread(fid, 32, '*char')');
    headers.slowchans(ii,1).channel = fread(fid, 1, 'int32');
    headers.slowchans(ii,1).ADFrequency = fread(fid, 1, 'int32');
    headers.slowchans(ii,1).adgain = fread(fid, 1, 'int32');
    headers.slowchans(ii,1).enabled = fread(fid, 1, 'int32');
    headers.slowchans(ii,1).preampgain = fread(fid, 1, 'int32');
    
    if(headers.version>=104)
        headers.slowchans(ii,1).spikechannel = fread(fid, 1, 'int32');
    else
        headers.slowchans(ii,1).spikechannel = 0;
        fread(fid, 1, 'int32');
    end

    if(headers.version>=105)
        headers.slowchans(ii,1).comment = deblank(fread(fid, 128, '*char')');
    else
        headers.slowchans(ii,1).comment = '';
        fread(fid, 128, 'char');
    end
    
    if(headers.version>=106)
        headers.slowchans(ii,1).srcid = fread(fid, 1, 'uchar');
        fread(fid, 1, 'uchar');
        headers.slowchans(ii,1).chanid = fread(fid, 1, 'ushort');
        fread(fid, 27, 'int32');
    else fread(fid, 28, 'int32');
    end
end

% Current position = 7504 + 1020*nchans + 296*nevchans + 296*nslowchans

% Keep track of the current position.
headers.datastart = ftell(fid);

% Unless instructed to read the entire file, stop here.
if(~fullread); return; end
if(DEBUG); fprintf(1,'Starting full file read.\n'); end

% If we are still going, then we were instructed to read the entire file to
% get accurate counts for spikes, events, and waveforms.

% Lets start by looking up the channel numbers.
if(headers.numDSPChannels>0)
    channums = vertcat(headers.chans.channel);
    maxchan = max(channums);
else
    maxchan = 0;
end

if(headers.numEventChannels>0)
    evnums = vertcat(headers.evchans.num);
    maxev = max(evnums);
else
    maxev = 0;
end

if(headers.numSlowChannels>0)
    slownums = vertcat(headers.slowchans.channel);
    maxslow = max(slownums);
else
    maxslow = 0;
end

tscounts   = zeros(27, maxchan);
wfcounts   = zeros(27, maxchan);
evcounts   = zeros(1,  maxev);
slowcounts = zeros(1,  maxslow+1);
slowfrags  = zeros(1,  maxslow+1);

% Read data headers and count number of each type of channel.
fseek(fid,headers.datastart,'bof');

numpointswave = headers.numPointsWave;

    function fillbuffer()
        dbufloc(1) = dbufloc-ndbuf;
        ndbuf(1) = min(floor((fsize-ftell(fid))/2),dbufsize);
        dbuf(1:ndbuf) = fread(fid,ndbuf,'*short');
    end

datablocks = 0;

dbufsize = 1024*1024*2;
dbuf = int16(zeros(dbufsize,1));
ndbuf = 0;
dbufloc = 1;

fillbuffer();

% type    = db(1);   db(end-7);
% upperts = db(2);   db(end-6);
% lowerts = db(3:4); db(end-5:end-4);
% chan    = db(5);   db(end-3);
% unit    = db(6);   db(end-2);
% nwaves  = db(7);   db(end-1);
% nwords  = db(8);   db(end);
dbh = int16(zeros(8,1));
needvals = 0;
havevals = 0;

while(ndbuf > 0)
    datablocks = datablocks + 1;
    if(DEBUG && mod(datablocks,100000)==0);
        fprintf('Block: %.0f\n',datablocks);
    end
    
    if(ndbuf-dbufloc>=7)
        dbh(1:8) = dbuf(dbufloc:dbufloc+7);
        dbufloc = dbufloc+8;
    else
        havevals(1) = max(0,ndbuf-dbufloc+1);
        dbh(1:havevals) = dbuf(dbufloc:ndbuf);
        needvals(1) = 8-havevals;
        dbufloc = dbufloc+havevals;
        fillbuffer();
        if(ndbuf-dbufloc+1<needvals); break; end
        dbh(havevals+1:8) = dbuf(dbufloc:dbufloc+needvals-1);
        dbufloc = dbufloc+needvals;
    end
    dbufloc = dbufloc+double(dbh(7)*dbh(8));

    if(dbh(1)==1)
        tscounts(dbh(6)+1,dbh(5)) = tscounts(dbh(6)+1,dbh(5)) + 1;
        wfcounts(dbh(6)+1,dbh(5)) = wfcounts(dbh(6)+1,dbh(5)) + dbh(7);
        if(dbh(8) ~= numpointswave(end))
            numpointswave(end+1) = dbh(8); %#ok<AGROW>
        end
    elseif(dbh(1) == 4)
        evcounts(dbh(5)) = evcounts(dbh(5)) + 1;
    elseif(dbh(1) == 5)
        slowfrags(dbh(5)+1) = slowfrags(dbh(5)+1) + 1;
        slowcounts(dbh(5)+1) = slowcounts(dbh(5)+1) + double(dbh(7)*dbh(8));
    end
end

if(~isscalar(numpointswave))
    headers.numPointsWave = max(numpointswave);
    
    warning('readPLXHeaders:NumPointsWaveMismatch',...
        'Number of points per wave is inconsistent or doesn''t match header information.');
end

headers.tscounts = tscounts;
headers.wfcounts = wfcounts;
headers.evcounts = evcounts;
headers.slowfrags = slowfrags;
headers.slowcounts = slowcounts;

for ii = 1:headers.numDSPChannels
    headers.chans(ii,1).nunits = find(tscounts(:,ii)>0,1,'last')-1;
end

end