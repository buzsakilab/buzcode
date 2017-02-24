function info = GetCurrentSession(verbose)

%GetCurrentSession - Get information about current session.
%
%  USAGE
%
%    GetCurrentSession
%    GetCurrentSession('verbose')  % detailed output
%
%    info = GetCurrentSession  % return info (do not display)

% Copyright (C) 2008-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA;
info = [];

verbose = nargin >= 1 && strcmp(lower(verbose),'verbose');

% Initialization
if isempty(DATA),
	disp('No session loaded');
	return
end

if nargout > 0,
	info.basename = DATA.session.basename;
	info.path = DATA.session.path;
	info.nChannels = DATA.nChannels;
	info.wideband = DATA.rates.wideband;
	info.lfp = DATA.rates.lfp;
	return
end

separator;
disp('Session');
disp(['  Name                           ' DATA.session.basename]);
disp(['  Disk Path                      ' DATA.session.path]);

separator;
disp('Channels');
disp(['  Number of Channels             ' int2str(DATA.nChannels)]);
disp(['  Wide-Band Sampling Rate        ' num2str(DATA.rates.wideband)]);
disp(['  LFP Sampling Rate              ' num2str(DATA.rates.lfp)]);

separator;
disp('Events');
disp(['  Total Number of Events         ' int2str(size(DATA.events.time,1))]);
ev = unique(DATA.events.description);
n = size(ev,1);
disp(['  Number of Different Types      ' int2str(n)]);
if verbose,
	for i = 1:n,
		if i < 10, white = ' '; else white = ''; end
		disp([' ' white int2str(i) '. Description                 ''' ev{i} '''']);
		n1 = length(GetEvents({ev{i}}));
		disp(['     Number                      ' int2str(n1)]);
	end
end
if n ~= 0,
	disp(['  Time of First Event            ' num2str(DATA.events.time(1))]);
	disp(['  Time of Last Event             ' num2str(DATA.events.time(end))]);
end

separator;
disp('Positions');
disp(['  Sampling Rate                  ' num2str(DATA.rates.video)]);
disp(['  Image Size (pixels)            ' int2str(DATA.maxX) 'x' int2str(DATA.maxY)]);
n = size(DATA.positions,1);
disp(['  Number of Position Samples     ' int2str(n)]);
if n ~= 0,
	disp(['  Number of Undetected Samples   ' int2str(sum(DATA.positions(:,1)==-1))]);
end

separator;
disp('Spikes');
nGroups = DATA.spikeGroups.nGroups;
disp(['  Number of Spike Groups         ' int2str(nGroups)]);
if verbose,
	for i = 1:nGroups,
		if i < 10, white1 = ' '; else white1 = ''; end
		clusters = DATA.spikes(DATA.spikes(:,2)==i,3);
		nClusters = length(unique(clusters));
		for j = 1:nClusters,
			if j < 10, white2 = ' '; else white2 = ''; end
			if j == 1, n = [white1 int2str(i) '.']; else n = '   ';end
			disp([' ' n ' Spikes in cluster ' white2 int2str(j) '        ' int2str(sum(clusters==j))]);
		end
	end
end

separator;

function separator

%  disp('-------------------------------------------------------------------------------------------------------');
disp(' ');