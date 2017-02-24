function SaveRippleEvents(filename,ripples,channelID)

%SaveRippleEvents - Save hippocampal ripple (~200Hz oscillations) events.
%
%  USAGE
%
%    SaveRippleEvents(filename,ripples,channelID)
%
%    filename       file to save to
%    ripples        ripple info as provided by <a href="matlab:help FindRipples">FindRipples</a>
%    channelID      channel ID (appended to the event description)
%
%  SEE
%
%    See also FindRipples, RippleStats, PlotRippleStats, SaveEvents.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 3,
  error('Incorrect number of parameters (type ''help <a href="matlab:help SaveRippleEvents">SaveRippleEvents</a>'' for details).');
end

n = size(ripples,1);
r = ripples(:,1:3)';
events.time = r(:);
for i = 1:3:3*n,
	events.description{i,1} = ['Ripple start ' int2str(channelID)];
	events.description{i+1,1} = ['Ripple peak ' int2str(channelID)];
	events.description{i+2,1} = ['Ripple stop ' int2str(channelID)];
end

SaveEvents(filename,events);