
function headers = checkPLXData(data, headers)
% checkPLXData - Check that PLX data and update/create PLX file headers.
%
% headers = checkPLXData(data, headers)
%
% checkPLXData is designed as a companion function to writePLXFile and
% readPLXHeaders. It takes spiking data that is destined for a PLX file,
% and (optionally) a header produced by readPLXHeaders, and checks that the
% data and header are ready for writing. This function is called by
% writePLXFile, so there is no need to call it first (although this won't
% harm anything). If headers is not provided, common reasonable defaults
% are used instead.
%
% WARNING: Continuous (slow channel) and event data is not yet supported by
% this function. All related data is thrown away and all related fields are
% zeroed out.
%
% Input:
%   data - A structure containing the spiking data
%     data.ts    - [n x 1] array of timestamps (seconds)
%     data.chans - [n x 1] array of channel assignments (1 based)
%     data.units - [n x 1] array of unit assignments (0 = unsorted)
%     data.waves - [n x npw] matrix of waveforms
%   headers - PLX file header information (optional). See readPLXHeaders.
%
% Output:
%   headers - PLX file header information, incorporating information from
%             data structure about spike and unit counts.
%
% See also writePLXFile, readPLXHeaders.
%
% AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% Copyright (c) 2011, Benjamin Kraus
% $Id$

idstr = 'BJK MATLAB R';
matlabrel = version('-release');

if(length(matlabrel)<=(18-length(idstr)))
    idstr = [idstr matlabrel];
elseif(length(matlabrel)>=15)
    idstr = ['BJK' matlabrel(1:15)];
else
    idstr = [idstr(1:18-length(matlabrel)) matlabrel];
end

% First verify that all fields within the data structure have the same
% number of entries.
if(size(data.ts,1) ~= size(data.chans,1) ||...
        size(data.ts,1) ~= size(data.units,1) ||...
        size(data.ts,1) ~= size(data.waves,1))
    error('All fields in PLX data must have same number of rows.');
end

headers.magic = 'PLEX'; % Required at the head of every PLX file.
headers.version = 106;  % This code is designed to write version 106.
headers.processingsoftware =  idstr;

headers.numPointsWave = size(data.waves,2);


% If headers are not provided, make them up.
if(nargin < 2)
    headers.ADFrequency = 40000;
    headers.lasttimestamp = round(max(data.ts)*headers.ADFrequency);
    headers.numDSPChannels = max(data.chans);
    headers.numPointsPreThr = 0;    
    headers.date = now();
    headers.fastread = 0;
    headers.waveformfreq = 40000;
    headers.trodalness = 4;             % Default for version < 103
    headers.datatrodalness = 1;         % Default for version < 103
    headers.bitsperspikesample = 12;    % Default for version < 103
    headers.bitsperslowsample = 12;     % Default for version < 103
    headers.spikeMaxMagnitudeMV = 3000; % Default for version < 103
    headers.slowMaxMagnitudeMV = 5000;  % Default for version < 103
    headers.spikePreAmpGain = 1000;     % Default for version < 105
    headers.acquiringsoftware = '';     % Default for version < 106
    headers.comment='';
    % Channel header information common all channels.
    chanheader.name = '';
    chanheader.signame = '';
    chanheader.channel = 0;
    chanheader.wfrate = 0;
    chanheader.sig = 0;
    chanheader.ref = 0;
    chanheader.gain = 1;
    chanheader.filter = 0;
    chanheader.threshold = 0;
    chanheader.method = 1;
    chanheader.nunits = 0;
    chanheader.template = zeros(5,64);
    chanheader.fit = zeros(5,1);
    chanheader.sortwidth = headers.numPointsWave;
    chanheader.boxes = zeros(5,2,4);
    chanheader.sortbeg = 0;
    chanheader.comment = ''; % Default for version < 105
    chanheader.srcid = '';   % Default for version < 106
    chanheader.chanid = 0;   % Default for version < 106

    % Create channel specific header values.
    for ch = 1:headers.numDSPChannels;
        DSPname = sprintf('sig%03u',ch);
        headers.chans(ch) = chanheader;
        headers.chans(ch).name = DSPname;
        headers.chans(ch).signame = DSPname;
        headers.chans(ch).channel = ch;
        headers.chans(ch).sig = ch;
    end
end

if isfield(headers,'numEventChannels') && isfield(headers,'numSlowChannels')&& ((headers.numEventChannels ~= 0 || headers.numSlowChannels ~= 0))
    wstr = ['Continuous (slow channel) and event data is not yet ',...
        'supported by this function. All related data is thrown away ',...
        'and all related fields are zeroed out.'];
    warning('checkPLXData:unsupported',wstr);
end
headers.numEventChannels = 0; % Event channels not supported
headers.numSlowChannels = 0;  % Slow channels not supported
headers.evcounts = zeros(1, 512);

% Make sure that numDSPChannels matches the number of channel headers.
headers.numDSPChannels = length(headers.chans);

% Shorthand for variables that will be used later.
nDSP = headers.numDSPChannels;
trodalness = headers.trodalness;

switch trodalness
    case 1; elstr = 'electrode';
    case 2; elstr = 'stereotrode';
    case 4; elstr = 'tetrode';
    otherwise; error('Unrecognized trodalness.');
end

% Initialize the timestamp, event, and waveform counts.
% tscounts and wvcounts must be at least [5 x 130]
tscounts = zeros(max(4,max(data.units))+1, max(nDSP,128)+2);

if(max(data.chans)>nDSP)
    error('More spike data channels than header spike channels.');
end

if(mod(nDSP,trodalness)~=0)
    error('Number of spike channel headers must be multiple of trodalness');
end

% Count the number of units on each channel and spikes in each unit.
nspks = zeros(nDSP,1);
tunits = [];
tref = 0;
for ch = 1:nDSP;
    tt = ceil(ch/trodalness);
    units = data.units(data.chans == ch);
    if(isempty(units))
        nspks(ch) = 0;
        headers.chans(ch).nunits = 0;
        tscounts(:,ch+1) = 0;
    else
        % Check if we have a unit assignment for this trode yet. If not,
        % save this unit assignment. If so, check that this is the same
        % unit assignment.
        if(tref ~= tt); tunits = units; tref = tt;
        elseif(length(tunits)~=length(units))
            error('All channels within a %s must have the same number of spikes.',elstr);
        elseif(~all(tunits == units))
            error('All channels within a %s must have the same unit assignments.',elstr);
        end
        nspks(ch) = length(units);
        headers.chans(ch).nunits = max(units);
        for jj = 0:max(data.units);
            tscounts(jj+1,ch+1) = sum(units == jj);
        end
    end
end
headers.tscounts = tscounts;
headers.wfcounts = tscounts;

% Check that each channel on a stereotrode or tetrode has the same number
% of spikes.
if(trodalness >= 2)
    nspks = reshape(nspks,trodalness,[])';
    maxspks = repmat(max(nspks,[],2),1,trodalness);
    nspks(nspks==0) = maxspks(nspks==0);
    if(any(min(nspks,[],2)~=max(nspks,[],2)))
        error('All channels within a %s must have the same number of spikes.',elstr);
    end
end