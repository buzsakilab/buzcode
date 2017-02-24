function spikes = GetSpikes(units,varargin)

%GetSpikes - Get spike timestamps.
%
%  NOTE
%
%    This function is provided for convenience. It simply calls <a href="matlab:help GetSpikeTimes">GetSpikeTimes</a>
%    using the same parameters. See this function for details.


% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

spikes = GetSpikeTimes(units,varargin{:});
