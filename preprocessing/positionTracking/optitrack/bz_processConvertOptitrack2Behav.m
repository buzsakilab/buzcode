function [tracking] = bz_processConvertOptitrack2Behav(fbasename,varargin)
% USAGE
%
%   bz_processConvertOptitrack2behavior(fbasename,varargin)
%    
% 
% INPUTS
%    fbasename   -basename of the recording (this is only used in generating
%                 the .mat output file)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'syncDatFile' name of binary file where sync signal is stored
%                   (default = filebasename_digitalin.dat)
%     'syncSampFq'  sampling freqeuncy of the sync signal (default= 20kHz)
%     'syncChan'    sync channel to read (default = 1)
%     'syncNbCh'    number of channels in the sync file (default = 1)
%     'posSampFq'   sampling frequency of the final pos file (after
%                   interpolation)
%     'columnOrder' order of output variales in .csv file: frame, time, rx,
%                   ry, rz, rw, x, y, z, error. Vector with postion of each
%                   variable in the order specified here
%    =========================================================================
%
% OUTPUTS
%
% 
% This function assumes that a Motive (Optitrack) behaviorioral session 
% has been exported to a .CSV file with default settings.
% (http://wiki.optitrack.com/index.php?title=Data_Export:_CSV)
%
% 
% David Tingley, 2017 (adapted from A. Peyrache, convert2pos.m)
% Antonio FR, 11/2018

syncDatFile = 'digitalin.dat'; % defaults if args aren't given... assuming a single digitalin channel on Intan
syncSampFq  = 20000;
syncChan    = 1;
syncNbCh    = 1;
posSampFq   = 120; %in Hz
columnOrder = 1:10; 

% Parse options
for i = 1:2:length(varargin)
  if ~isa(varargin{i},'char')
    error(['Parameter ' num2str(i+3) ' is not a property .']);
  end
  switch(lower(varargin{i}))
    case 'syncdatfile'
      syncDatFile = varargin{i+1};
      if ~isa(duration,'char')
        error('Incorrect value for property ''syncDatFile'' .');
      end
    case 'syncsampfq'
      syncSampFq = varargin{i+1};
      if ~isa(syncSampFq,'numeric')
        error('Incorrect value for property ''syncsampfq'' .');
      end
    case 'syncchan'
      syncChan = varargin{i+1};
      if ~isa(start,'numeric')pw
        error('Incorrect value for property ''syncChan'' .');
      end
        if start < 0, start = 0; end
    case 'syncnbch'
      syncNbCh = varargin{i+1};
      if ~isa(syncNbCh,'numeric')
        error('Incorrect value for property ''syncNbCh'' .');
      end
    case 'possampfq'
      posSampFq = varargin{i+1};
      if ~isa(posSampFq,'numeric')
        error('Incorrect value for property ''posSampFq'' .');
      end
    case 'columnorder'
      columnOrder = varargin{i+1};
      if ~isa(columnOrder,'numeric')
        error('Incorrect value for property ''syncNbCh'' .');
      end      
    otherwise
      error(['Unknown property ''' num2str(varargin{i}) ''' .']);
  end
end


if exist([fbasename '.csv'],'file') % assumes csv naming matches basename
    dat = importdata([fbasename '.csv']);
    dat = scrubTracking(dat); % smooth and interpolate tracking data
elseif ~isempty(dir('*.csv')) % looks for any csv file in the current recording folder
    csv = dir('*.csv');
    dat = importdata(csv.name);
    dat = scrubTracking(dat);
elseif ~isempty(dir('Session*')) % Looks for a /Session*/ folder that Motive/Optitrack has created
    d = dir('Session*');
    cd(d.name)
    csv = dir('*.csv');
    dat = importdata(csv.name);
    dat = scrubTracking(dat);
    cd ..
else
    error('couldnt find a .csv tracking file')
end

%%
pos = dat.data; % all tracking variables   

    % read optitrack sync channel
    fid = fopen(syncDatFile); 
    dig = fread(fid,[syncNbCh inf],'int16=>int16');  % default type for Intan digitalin
    dig = dig(syncChan,:);
    t = (0:length(dig)-1)'/syncSampFq; % time vector in sec
    
if exist('digitalIn.events.mat','file') 
   % new way of getting frames, better than before. Requiere running before
   % bz_getDigitalIn
   load('digitalIn.events.mat');
   temp = digitalIn.timestampsOn{syncChan}; % frame onsets
   frameT = temp+(0.5/posSampFq); % center of the frame
   
else
    % this is the old way. Prone to errors when there are more that digital input. 
    dPos = find(diff(dig)==1);
    dNeg = find(diff(dig)==-1);

    if length(dPos) == length(dNeg)+1
        dPos = dPos(1:end-1);
    elseif length(dNeg) == length(dPos)+1
        dNeg = dNeg(1:end-1);
    elseif abs(length(dNeg)-length(dPos)) > 1
        warning('some problem with frames');
        keyboard
    end
    % Frame timing is the middle of shuter opening
    frameT = (t(dPos)+t(dNeg))/2;
end

% The system sometimes (rarely) keeps on recording a few frames after software stopped
% recording. So we skip the last frames of the TTL

if length(frameT)<size(pos,1)
    warning('Too many video frames!'); % maybe because intan was stopped before Optitrack
    %keyboard
    pos(pos==-1) = NaN;
    pos = pos(1:length(frameT),:); % ???
elseif length(frameT)>size(pos,1)
    frameT = frameT(1:size(pos,1));
end

% We now interpolate the data at 120 Hz (or sampling fqcy specified in arguments)
recDuration = length(dig)/syncSampFq;

timestamps = (0:1/posSampFq:recDuration-1/posSampFq)';
pos(pos==-1) = NaN;
try 
newPos = interp1(frameT,pos,timestamps);
catch
end
%error('need to remove NaN garbage still')

%% create output structure 

tracking.timestamps = timestamps;
tracking.frameCount = newPos(:,columnOrder(2));

tracking.orientation.rx = newPos(:,columnOrder(3));
tracking.orientation.ry = newPos(:,columnOrder(4));
tracking.orientation.rz = newPos(:,columnOrder(5));
tracking.orientation.rw = newPos(:,columnOrder(6));

tracking.position.x = newPos(:,columnOrder(7));
tracking.position.y = newPos(:,columnOrder(8));
tracking.position.z = newPos(:,columnOrder(9));

if  size(newPos,2) >= 10
tracking.errorPerMarker = newPos(:,columnOrder(10));
end

save([fbasename '.tracking.behavior.mat'],'tracking');

end

