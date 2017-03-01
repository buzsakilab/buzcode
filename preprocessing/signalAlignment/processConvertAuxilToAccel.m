function [behav] = processConvertAuxilToAccel(fbasename,varargin) 
% USAGE
%
%   processConvertAuxilToAccel(fbasename,varargin)
%    
% 
%  INPUTS
%    fbasename   basename of the dat file
%    
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'syncDatFile' name of binary file where sync signal is stored
%                   (default = filebasename_digitalin.dat)
%     'accelFile'   basename of .dat with accelerometry data
%     'recSampFq'  sampling freqeuncy of the sync signal (default= 20kHz)
%     'syncChan'    sync channel to read (default = 1)
%     'numChans'    number of channels in the sync file (default = 1)
%     'accelSampFq'   sampling frequency of the final pos file (after
%                   interpolation)
%    =========================================================================
%
% OUTPUT
%    behav a behavior struct with the following fields added
%    behav.accel.x - acceleration in the x direction
%    behav.accel.y - acceleration in the x direction
%    behav.accel.z - acceleration in the x direction
%    behav.accel.ts - timestamps of downsamples acceleration signal


accelFile = 'auxiliary.dat'; % default intan file
recSampFq  = 20000;
numChans = 3;
accelSampFq = 1250;

for i = 1:2:length(varargin)
  if ~isa(varargin{i},'char')
    error(['Parameter ' num2str(i+3) ' is not a property .']);
  end
  switch(lower(varargin{i}))
    case 'accelfile'
      syncDatFile = varargin{i+1};
      if ~isa(duration,'char')
        error('Incorrect value for property ''syncDatFile'' .');
      end
    case 'recsampFq'
      recSampFq = varargin{i+1};
      if ~isa(frequency,'numeric')
        error('Incorrect value for property ''recSampFq'' .');
      end
    case 'numchans'
      numChans = varargin{i+1};
      if ~isa(numChans,'numeric')
        error('Incorrect value for property ''numChans'' .');
      end
    case 'accelsampfq'
      posSampFq = varargin{i+1};
      if ~isa(posSampFq,'numeric')
        error('Incorrect value for property ''posSampFq'' .');
      end
    otherwise
      error(['Unknown property ''' num2str(varargin{i}) ''' .']);
  end
end

if ~exist(accelFile)
    error(['couldnt find ' accelFile ' file for accelerometer data...']) 
end
    % read data into workspace
    fid = fopen(accelFile); 
    dat = fread(fid,[numChans inf],'uint16=>uint16');
    fclose(fid);

%     fid = fopen('time.dat');
%     ts = fread(fid,[1,inf],'int32=>int32');  % we can do this without IO
    ts = 0:1/recSampFq:length(dat)-1/recSampFq;
%     fclose(fid);

    % downsample to 1250 Hz
    behav.accel.x = uint16(resample((double(dat(1,:))),accelSampFq,recSampFq));
    behav.accel.y = uint16(resample((double(dat(2,:))),accelSampFq,recSampFq));
    behav.accel.z = uint16(resample((double(dat(3,:))),accelSampFq,recSampFq));
    behav.accel.ts = ts;
end