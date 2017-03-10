function stats = MapStats(map,varargin)

%MapStats - Compute statistics for a map Z = f(X,Y), or a curve Z = f(X).
%
%  Compute statistics on a continuous map, where one time-varying variable Z
%  is represented as a function of one or two time-varying variables X and Y.
%  The variable Z can either be a point process (typically, a list of spike
%  timestamps) or a continuous measure (e.g. the instantaneous velocity of
%  the animal, the spectral power of an LFP channel in a given frequency band,
%  the coherence between two oscillating LFP channels, etc.) Typical examples
%  of X and Y include spatial coordinates and angular directions.
%
%  USAGE
%
%    stats = MapStats(map,<options>)
%
%    map            map obtained using <a href="matlab:help Map">Map</a>
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'threshold'   values above threshold*peak belong to the field
%                   (default = 0.2)
%     'minSize'     fields smaller than this size are considered spurious
%                   and ignored (default = 100 for 2D, 10 for 1D)
%     'minPeak'     peaks smaller than this size are considered spurious
%                   and ignored (default = 1)
%     'type'        'll' if X and Y are linear, 'cl' if X is circular and Y
%                   linear, 'lc' if X is linear and Y circular, or 'cc' if X
%                   and Y are circular - for 1D data, a single letter is used
%                   (default = 'll')
%     'verbose'     display processing information (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    stats.x             abscissa of the maximum value (in bins)
%    stats.y             ordinate of the maximum value (in bins)
%    stats.peak          in-field maximum value
%    stats.mean          in-field mean value
%    stats.size          field size (in bins)
%    stats.field         field (1 = bin in field, 0 = bin not in field)
%    stats.fieldX        field x boundaries (in bins)
%    stats.fieldY        field y boundaries (in bins)
%    stats.specificity   spatial specificity (Skaggs et al., 1993)
%
%    For 1D circular data:
%
%    stats.m             mean angle
%    stats.mode          distribution mode (in bins)
%    stats.r             mean resultant length
%    stats.k             von Mises concentration
%
%  SEE
%
%    See also Map, FiringMap, FiringCurve.

% Copyright (C) 2002-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).');
end

if ~isfield(map,'z'),
	map.z = map.rate;
	map = rmfield(map,'rate');
end

% Default values
threshold = 0.2;
minPeak = 1;
type = 'll';
verbose = 0;

nDims = sum(size(map.z)>=2);
if nDims == 2,
	minSize = 100;
else
	minSize = 10;
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).']);
	end
	switch(lower(varargin{i})),

		case 'threshold',
			threshold = varargin{i+1};
			if ~isdscalar(threshold,'>=0','<=1'),
				error('Incorrect value for property ''threshold'' (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).');
			end

		case 'minsize',
			minSize = varargin{i+1};
			if ~isdscalar(minSize,'>=0'),
				error('Incorrect value for property ''minSize'' (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).');
			end

		case 'minpeak',
			minPeak = varargin{i+1};
			if ~isdscalar(minPeak,'>=0'),
				error('Incorrect value for property ''minPeak'' (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).');
			end

		case 'type',
			type = lower(varargin{i+1});
			if (nDims == 1 && ~isstring_FMAT(type,'c','l')) || (nDims == 2 && ~isstring_FMAT(type,'cc','cl','lc','ll')),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).');
			end

		case {'verbose','debug'},
			verbose = lower(varargin{i+1});
			if ~isstring_FMAT(verbose,'on','off'),
				error('Incorrect value for property ''verbose'' (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).');
			end
			verbose = strcmp(verbose,'on');

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).']);

  end
end

% Are X and/or Y circular?
circX = size(map.z,2) > 1 && strcmp(type(1),'c');
circY = size(map.z,1) > 1 && ((size(map.z,2) > 1 && strcmp(type(2),'c')) || strcmp(type(1),'c'));

% Default values
stats.x = NaN;
stats.y = NaN;
stats.field = logical(zeros(0,0,0));
stats.size = 0;
stats.peak = 0;
stats.mean = 0;
stats.fieldX = [NaN NaN];
stats.fieldY = [NaN NaN];
stats.specificity = 0;
stats.m = nan;
stats.r = nan;
stats.mode = nan;
stats.k = nan;

if isempty(map.z), return; end

% Compute the spatial specificity of the map, based on the formula proposed by Skaggs et al. (1993).
T = sum(sum(map.time));
if T == 0,
  stats.specificity = 0;
else
  occupancy = map.time/(T+eps);
  m = sum(sum(map.count))/(sum(sum(map.time)+eps));
  if m == 0,
    stats.specificity = 0;
  else
    logArg = map.count/m;
    logArg(logArg <= 1) = 1;

%      stats.specificity = sum(sum(map.count.smooth.*log2(logArg).*occupancy))/m;
    stats.specificity = sum(sum(map.count.*log2(logArg).*occupancy))/m;
  end
end

% Determine the field as the connex area around the peak where the value or rate is > threshold*peak
% There can be two or more fields

if max(max(map.z)) == 0,
  stats.field = logical(zeros(size(map.z)));
  return;
end

nBinsX = max([1 length(map.x)]);	% minimum number of bins is 1
nBinsY = max([1 length(map.y)]);

% Each time we find a field, we will remove it from the map; make a copy first
z = map.z;
% Try to find more fields until no remaining bin exceeds min value
i = 1;
while true,
	% Are there any candidate (unvisited) peaks left?
	[peak,idx] = max(z(:));
	if peak < minPeak, break; end
	% Determine coordinates of largest candidate peak
	[y,x] = ind2sub(size(z),idx);
	% Find field (using min threshold for inclusion)
	field1 = FindField(z,x,y,peak*threshold,circX,circY);
	size1 = sum(field1(:));
	% Does this field include two coalescent subfields?
	% To answer this question, we simply re-run the same field-searching procedure on the field
	% we then either keep the original field or choose the subfield if the latter is less than
	% 1/2 the size of the former
	m = peak*threshold;
	field2 = FindField(z-m,x,y,(peak-m)*threshold,circX,circY);
	size2 = sum(field2(:));
	if size2< 1/2*size1,
		field = field2;
		tc = ' ';sc = '*'; % for debugging messages
	else
		field = field1;
		tc = '*';sc = ' '; % for debugging messages
	end
	% Display debugging info
	if verbose,
		disp([int2zstr(i,2) ') peak  ' num2str(peak) ' @ (' int2str(x) ',' int2str(y) ')']);
		disp([' ' tc ' field size       ' int2str(size1)]);
		disp([' ' sc ' subfield size    ' int2str(size2)]);
		disp(' ');
		figure;
		if nDims == 1,
			plot(z);hold on;
			PlotIntervals(ToIntervals(field1),'rectangles');
			PlotIntervals(ToIntervals(field2),'bars');
			ylabel(tc);
		else
			subplot(3,1,1);imagesc(z);xlabel('Data');
			subplot(3,1,2);imagesc(field1);clim([0 1]);xlabel('Field');
			subplot(3,1,3);imagesc(field2);clim([0 1]);ylabel(tc);xlabel('Subfield');
		end
	end
	fieldSize = sum(field(:));
	% Keep this field if its size is sufficient
	if fieldSize > minSize,
		stats.field(:,:,i) = field;
		stats.size(i) = fieldSize;
		stats.peak(i) = peak;
		stats.mean(i) = mean(z(field));
		idx = find(field & z == peak);
		[stats.y(i),stats.x(i)] = ind2sub(size(z),idx(1));
		[x,y] = FieldBoundaries(field,circX,circY);
		[stats.fieldX(i,:),stats.fieldY(i,:)] = FieldBoundaries(field,circX,circY);
		i = i + 1;
	end
	% Mark field bins as visited
	z(field) = NaN;
	if all(isnan(z)), break; end
end

% ----------------------------- Circular statistics -----------------------------

if ~(nDims == 1 && circX), return; end

complex = map.z .* exp(j*map.x*2*pi) / sum(map.z);
stats.m = angle(nanmean(complex));
stats.r = abs(nansum(complex));
%  [~,stats.mode] = max(map.z);
[maxZ, iMax] = max(map.z);
TMax = map.x(iMax);
stats.mode = TMax*2*pi;

% Compute kappa (Fisher, "Statistical Analysis of Circular Data" p.88)
if stats.r < 0.53,
	stats.k = 2*stats.r+stats.r^3+5*stats.r^5/6;
elseif stats.r < 0.85,
	stats.k = -0.4+1.39*stats.r+0.43/(1-stats.r);
else
	stats.k = 1/(stats.r^3-4*stats.r^2+3*stats.r);
end

% ------------------------------- Helper functions -------------------------------

% Field boundaries (circumscribed rectangle)

function [x,y] = FieldBoundaries(field,circX,circY)

% Find boundaries
x = find(any(field,1));
if isempty(x),
	x = [NaN NaN];
else
	x = [x(1) x(end)];
end
y = find(any(field,2));
if isempty(y),
	y = [NaN NaN];
else
	y = [y(1) y(end)];
end

% The above works in almost all cases; it fails however for circular coordinates if the field extends
% around an edge, e.g. for angles between 350° and 30°

if circX && x(1) == 1 && x(2) == size(field,2),
	xx = find(~all(field,1));
	if ~isempty(xx),
		x = [xx(end) xx(1)];
	end
end
if circY && y(1) == 1 && y(2) == size(field,1),
	yy = find(~all(field,2));
	if ~isempty(yy),
		y = [yy(end) yy(1)];
	end
end

