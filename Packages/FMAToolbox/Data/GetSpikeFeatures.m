function [data,indices] = GetSpikeFeatures(group)

%GetSpikeFeatures - Get spike features.
%
%  USAGE
%
%    features = GetSpikeFeatures(group)
%
%    group          electrode group
%
%  OUTPUT
%
%    features       (timestamp,group,cluster,features...) t-uples.
%
%  NOTE
%
%    Contrary to e.g. GetSpikes, this function retrieves information for all clusters
%    of a given electrode group, rather than from a list of selected units. This is
%    because .fet files are usually very large and it would be very inefficient to
%    repeatedly read the same file for each cluster. It is preferable to load the whole
%    file only once, and then select individual clusters if appropriate.
%
%  SEE
%
%    See also LoadSpikeFeatures, IsolationQuality.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA;
if isempty(DATA),
	error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help GetSpikeFeatures">GetSpikeFeatures</a>'' for details).');
end

if ~isiscalar(group,'>=0'),
	error('Incorrect unit (type ''help <a href="matlab:help GetSpikeFeatures">GetSpikeFeatures</a>'' for details).');
end

% Load all waveforms for this group
filename = [DATA.session.path '/' DATA.session.basename '.fet.' int2str(group)];
rate = DATA.rates.wideband;
data = LoadSpikeFeatures(filename,rate);
