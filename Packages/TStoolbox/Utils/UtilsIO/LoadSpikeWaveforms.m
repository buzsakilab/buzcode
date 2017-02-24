%LoadSpikeWaveforms - Read spike waveforms from file.
%
%  USAGE
%
%    waveforms = LoadSpikeWaveforms(filename,nChannels,nSamples)
%
%    filename            spike waveform file name
%    nChannels           number of channels in electrode group
%    nSamples            number of samples per waveform

% Copyright (C) 2004-2006 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

function waveforms = LoadSpikeWaveforms(filename,nChannels,nSamples)

if nargin < 3,
  error('Incorrect number of parameters (type ''help LoadSpikeWaveforms'' for details).');
end
waveforms = LoadBinary(filename);
waveforms = reshape(waveforms,nChannels,nSamples,[]);