function Process_ConvertOptitrack2Pos(fbasename,varargin)

% USAGE
%
%   Process_ConvertOptitrack2Pos(fbasename,varargin)
%    
% 
% 
%    fbasename   basename of the video file (should be filebasenmae.avi)
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
% 
% This function assumes that a Motive (Optitrack) behavioral session 
% has been exported to a .CSV file with default settings.
% (http://wiki.optitrack.com/index.php?title=Data_Export:_CSV)
%
% 
%

if strcmp(fbasename,'')
    d = dir('*sessionInfo.mat');
    load(d.name);
    fbasename = sessionInfo.FileName;
end

warning('this functions is now deprecated and has been replaced by processConvertOptitrack2Behav.m')

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
    if length(csv) > 1
        dat=[];
       for i=1:length(csv)
           dat = [dat; importdata(csv(i).name)];
       end
    else
        dat = importdata(csv.name);
    end
    
    if isstruct(dat)
        d=[];
        for i=1:length(csv)
             d = [d;dat(i).data]; 
        end
        clear dat; dat.data = d;
    end
    dat = scrubTracking(dat);
elseif ~isempty(dir('Session*')) % Looks for a /Session*/ folder that Motive/Optitrack has created
    d = dir('Session*');
    cd(d.name)
    csv = dir('*.csv');
    if length(csv) > 1
        dat=[];
       for i=1:length(csv)
           dat = [dat; importdata(csv(i).name)];
       end
    else
        dat = importdata(csv.name);
    end
    if isstruct(dat)
        d=[];
        for i=1:length(csv)
             d = [d;dat(i).data]; 
        end
        clear dat; dat.data = d;
    else
        d.data = dat;
        dat = d;
    end
    dat = scrubTracking(dat);
    cd ..
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
% recDuration = length(dig)/syncSampFq;
% 
% timestamps = (0:1/posSampFq:recDuration-1/posSampFq)';
% pos(pos==-1) = NaN;
% newPos = interp1(frameT,pos,timestamps,'linear');

% timestamps(isnan(newPos(:,2))) = [];
% newPos(isnan(newPos(:,2)),:)=[];

dlmwrite([fbasename '.pos'],[frameT pos],'delimiter','\t', 'precision', 32);

end

