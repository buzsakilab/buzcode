function processConvertOptitrack2Behav(fbasename,varargin)
% USAGE
%
%   bz_processConvertOptitrack2Behav(fbasename,varargin)
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
%    =========================================================================
%
% OUTPUTS
%
% 
% This function assumes that a Motive (Optitrack) behavioral session 
% has been exported to a .CSV file with default settings.
% (http://wiki.optitrack.com/index.php?title=Data_Export:_CSV)
%
% 
% David Tingley, 2017 (adapted from A. Peyrache, convert2pos.m)

syncDatFile = ['digitalin.dat']; % defaults if args aren't given... assuming a single digitalin channel on Intan
syncSampFq  = 20000;
syncChan    = 1;
syncNbCh    = 1;
posSampFq   = 120; %in Hz

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
      if ~isa(frequency,'numeric')
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
    otherwise
      error(['Unknown property ''' num2str(varargin{i}) ''' .']);
  end
end


if exist([fbasename '.csv'],'file') % assumes csv naming matches basename
    dat = importdata([fbasename '.csv']);
    dat = scrubTracking(dat); 
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

pos = dat.data;    

fid = fopen(syncDatFile); 
dig = fread(fid,[syncNbCh inf],'int16=>int16');  % default type for Intan digitalin
dig = dig(syncChan,:);

t = (0:length(dig)-1)'/syncSampFq;

dPos = find(diff(dig)==1);
dNeg = find(diff(dig)==-1);

if length(dPos) == length(dNeg)+1
    dPos = dPos(1:end-1);
elseif length(dPos) > length(dNeg)+1 || length(dPos) < length(dNeg)
     keyboard
end

% Frame timing is the middle of shuter opening
frameT = (t(dPos)+t(dNeg))/2;

% The system sometimes (rarely) keeps on recording a few frames after software stopped
% recording. So we skip the last frames of the TTL

if length(frameT)<size(pos,1)
    warning('Too many video frames!')
    keyboard
elseif length(frameT)>size(pos,1)
    frameT = frameT(1:size(pos,1));
end

% We now interpolate the data at 120 Hz (or sampling fqcy specified in
% arguments)
recDuration = length(dig)/syncSampFq;

timestamps = (0:1/posSampFq:recDuration-1/posSampFq)';
pos(pos==-1) = NaN;
newPos = interp1(frameT,pos,timestamps);

error('need to remove NaN garbage still')

behav.position.x = pos(:,2);
behav.position.y = pos(:,3);
behav.position.z = pos(:,4);
behav.orientation.rx = pos(:,5);
behav.orientation.ry = pos(:,6);
behav.orientation.rz = pos(:,7);
behav.orientation.rw = pos(:,8);
behav.timestamps = pos(:,1);
behav.errorPerMarker = pos(:,9);
behav.frameCount = pos(:,10);
save([fbasename '.tracking.behavior.mat'],'behav')

end

