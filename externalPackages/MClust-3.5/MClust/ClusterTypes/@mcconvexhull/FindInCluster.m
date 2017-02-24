function [f, MCC] = FindInCluster(MCC)

% f = FindInCluster(MCC, data)
%
% INPUTS
%     MCC - a MCCluster
%     data - nS x nD data
%
% OUTPUTS
%     f - indices of all points in the cluster
%
% ADR 1999
% version 1.2
%
% Status: PROMOTED (Release version)
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

% Modified ncst 02 March 03, removed xlbls and ylbls, added cowen's mod for
%   ForbiddenPoints, and to use the MClust_FDdn if the current directory is
%   changed.

nBounds = length(MCC.xdimNames);
if nBounds ~= length(MCC.ydimNames)
	error('Unaligned data in FindInCluster');
end
if nBounds ~= length(MCC.cx)
	error('Unaligned data in FindInCluster');
end
if nBounds ~= length(MCC.cy)
	error('Unaligned data in FindInCluster');
end

if nBounds == 0
	f = [];
else
	for iB = 1:nBounds

	   % X
	   temp = load(MCC.xdimSources{iB,1},'-mat');
	   data(:,1) = temp.FeatureData(:,MCC.xdimSources{iB,2});

	   % Y
	   temp = load(MCC.ydimSources{iB,1},'-mat');
	   data(:,2) = temp.FeatureData(:,MCC.ydimSources{iB,2});

		if exist('f', 'var')
			f0 = InPolygon(data(f,1), data(f,2), MCC.cx{iB}, MCC.cy{iB});
			f = f(logical(f0));
		else
			f = find(InPolygon(data(:,1), data(:,2), MCC.cx{iB}, MCC.cy{iB}));
		end
		if isempty(f)
			break;
		end
	end

	f = sort(f);
end