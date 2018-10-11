% ConvertBasler2pos - Convert Basler data to pos file format
%
%  USAGE
%
%    ConvertBasler2pos(filebasename,<options>)
%
%    filebasename   basename of the video file (should be filebasenmae.avi)
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
% DEPENDENCIES:
%
%   bz_LoadBinary
%   Process_DetectLED


% Copyright (C) 2015 Adrien Peyrache
% edited by David Tingley, 2017
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
function Process_ConvertBasler2Pos(fbasename,varargin)

warning('This function is now deprecated and has been replaced by processConvertLED2Behav.m')

syncDatFile = ['analogin.dat'];
syncSampFq  = 20000;
syncChan    = 1;
syncNbCh    = 1;
posSampFq   = 30; %in Hz

% Parse options
for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help Process_ConvertBasler2Pos'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'syncdatfile',
      syncDatFile = varargin{i+1};
      if ~isa(duration,'char')
        error('Incorrect value for property ''syncDatFile'' (type ''help Process_ConvertBasler2Pos'' for details).');
      end
    case 'syncsampfq',
      syncSampFq = varargin{i+1};
      if ~isa(frequency,'numeric')
        error('Incorrect value for property ''syncsampfq'' (type ''help Process_ConvertBasler2Pos'' for details).');
      end
    case 'syncchan',
      syncChan = varargin{i+1};
      if ~isa(start,'numeric')
        error('Incorrect value for property ''syncChan'' (type ''help Process_ConvertBasler2Pos'' for details).');
      end
		if start < 0, start = 0; end
    case 'syncnbch',
      syncNbCh = varargin{i+1};
      if ~isa(syncNbCh,'numeric')
        error('Incorrect value for property ''syncNbCh'' (type ''help Process_ConvertBasler2Pos'' for details).');
      end
    case 'possampfq',
      posSampFq = varargin{i+1};
      if ~isa(posSampFq,'numeric')
        error('Incorrect value for property ''posSampFq'' (type ''help Process_ConvertBasler2Pos'' for details).');
      end
    
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help Process_ConvertBasler2Pos'' for details).']);
  end
end




    if ~exist([fbasename '.led'])    
%         for thresh = [.05 .15 .3 .7 .9 1.5 5 10]
            thresh = .3;
        Process_DetectLED(fbasename,thresh);
%         end
    end
    pos = load([fbasename '.led']);
    
    
    fid = fopen(syncDatFile); 
    dat = fread(fid,[syncNbCh inf],'uint16=>uint16');
    dat = dat(syncChan,:);
    
    % the line below does not work correctly with Intan analogin.dat files
%     dat = bz_LoadBinary(syncDatFile,'nchannels',syncNbCh,'channels',syncChan);

    t = (0:length(dat)-1)'/syncSampFq;
    
    %It's a pull-off, trigger is toward 0
%     dat = double(fastrms(fastrms(fastrms(dat,10),50),100)<500);
	dat = uint16(smooth(single(dat)))<20000;

% if bad camera pulses, use below
%     dat = smooth(dat,100);
%     dat = dat<0;


    dPos = find(diff(dat)==1);
    dNeg = find(diff(dat)==-1);
    
    if length(dPos) == length(dNeg)+1
        dPos = dPos(1:end-1);
    elseif length(dPos) > length(dNeg)+1 || length(dPos) < length(dNeg)
         keyboard
    end
    
    %Frame timing is the middle of shuter opening
    frameT = (t(dPos)+t(dNeg))/2;
    
    %The camera keeps on recording a few frames after software stopped
    %recording. So we skipthe last frames of the TTL
    if length(frameT)<size(pos,1);
        warning('Too many video frames!')
        pause
    elseif length(frameT)>size(pos,1);
        frameT = frameT(1:size(pos,1));
    end
    
    %We now interpolate the data at 60 Hz (or sampling fqcy specified in
    %arguments)
    recDuration = length(dat)/syncSampFq;
    
    timestamps = (0:1/posSampFq:recDuration-1/posSampFq)';    
%     pos(pos==-1) = NaN;
    newPos = interp1(frameT,pos,timestamps);
    
    dlmwrite([fbasename '.pos'],[timestamps newPos],'delimiter','\t', 'precision', 32);

end

