%LoadBinary - Load data from a binary file.
%
%  USAGE
%
%    data = LoadBinary(filename,<options>)
%
%    filename       file to read
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'duration'    duration to read (in s) (default = Inf)
%     'frequency'   sampling rate (in Hz) (default = 20kHz)
%     'start'       position to start reading (in s) (default = 0)
%     'nChannels'   number of data channels in the file (default = 1)
%     'channels'    channels to read (default = all)
%     'precision'   sample precision (default = 'int16')
%     'skip'        number of bytes to skip after each value is read
%                   (default = 0)
%    =========================================================================

% Copyright (C) 2004-2006 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

function data = LoadSpikeWaveF(filename,nChannels,nSamples,CluIx)

% Default values
start = 0;
precision = 'int16';
skip = 0;
duration = Inf;
frequency = 20000;
channels = [];

if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end
f = fopen(filename,'r');
if f == -1,
	error(['Cannot read ' filename ' (insufficient access rights?).']);
end

nSamplesPerSpk = nChannels*nSamples;

data = [];
chunkSz = 10000;
nbChunks = floor(CluIx(end)/chunkSz);
data = zeros(nSamplesPerSpk*length(CluIx),1,'int16');
dIx=1;

for ii=1:nbChunks
    ix = CluIx>(ii-1)*chunkSz & CluIx<=ii*chunkSz;
    if length(ix)
    d = fread(f,nSamplesPerSpk*chunkSz,precision);
    d = reshape(d,nChannels,nSamples,[]);
    d = d(:,:,CluIx(ix)-chunkSz*(ii-1));
    d = d(:);
    data(dIx:dIx-1+sum(ix)*nSamplesPerSpk) = d;
    dIx = dIx+sum(ix)*nSamplesPerSpk;
    end
end

if any(CluIx>chunkSz*nbChunks)
    ix = CluIx>chunkSz*nbChunks;
    d = fread(f,nSamplesPerSpk*(CluIx(end)-chunkSz*nbChunks),precision);
    d = reshape(d,nChannels,nSamples,[]);
    d = d(:,:,CluIx(ix)-chunkSz*nbChunks);
    d = d(:);
    data(dIx:dIx-1+sum(ix)*nSamplesPerSpk) = d;
end

fclose(f);
data = reshape(data,nChannels,nSamples,[]);

