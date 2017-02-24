function AddSpikeTimes(spikes)

%AddSpikeTimes - Add spike timestamps.
%
%  USAGE
%
%    AddSpikeTimes(spikes)
%
%  SEE
%
%    See also ThresholdSpikes.


% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA;
if isempty(DATA),
	error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end

% Parse parameters
if ~issamples(spikes,'#2'),
	error('Incorrect spikes (type ''help <a href="matlab:help AddSpikeTimes">AddSpikeTimes</a>'' for details).');
end

DATA.spikes = [DATA.spikes;spikes];
