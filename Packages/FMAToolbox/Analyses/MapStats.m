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
%     'debug'       display processing information (default = 'off')
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
%  SEE
%
%    See also Map, FiringMap, FiringCurve.

% Copyright (C) 2002-2011 by Michaël Zugaro
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
debug = 0;

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

		case 'debug',
			debug = lower(varargin{i+1});
			if ~isstring_FMAT(debug,'on','off'),
				error('Incorrect value for property ''debug'' (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).');
			end
			debug = strcmp(debug,'on');

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
stats.field = logical([]);
stats.size = 0;
stats.peak = 0;
stats.mean = 0;
stats.fieldX = [NaN NaN];
stats.fieldY = [NaN NaN];

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
	% Are there any candidate (unvisited) bins left?
	candidates = z > minPeak;
	if ~any(candidates(:)), break; end
	% Determine coordinates of first candidate bin
	idx = find(candidates);
	[y1,x1] = ind2sub(size(candidates),idx(1));
	% Find field (at this point, all bins >0 are counted in the field)
	field1 = FindField(z,x1,y1,0,circX,circY);
	size1 = sum(field1(:));
	% Compute peak value and location and refine field using min threshold for inclusion
	peak = max(z(field1));
	idx = find(field1 & z == peak);
	[y2,x2] = ind2sub(size(z),idx(1));
	field2 = FindField(z,x2,y2,peak*threshold,circX,circY);
	size2 = sum(field2(:));
	% Does this field include two coalescent subfields?
	% To answer this question, we simply re-run the same field-searching procedure on the field
	% we then either keep the original field or choose the subfield if the latter is less than
	% 1/2 the size of the former
	m = peak*threshold;
	field3 = FindField(z-m,x2,y2,(peak-m)*threshold,circX,circY);
	size3 = sum(field3(:));
	if size3 < 1/2*size2,
		field = field3;
		tc = ' ';sc = '*'; % for debugging messages
	else
		field = field2;
		tc = '*';sc = ' '; % for debugging messages
	end
	% Display debugging info
	if debug,
		disp([int2str(i) ') candidates        ' int2str(sum(candidates(:)))]);
		disp(['   starting          (' int2str(x1) ',' int2str(y1) ')']);
		disp(['   raw field size    ' int2str(size1)]);
		disp([' ' tc ' thresholded size  ' int2str(size2)]);
		disp([' ' sc ' subfield size     ' int2str(size3)]);
		disp(['   peak              ' num2str(peak) ' @ (' int2str(x2) ',' int2str(y2) ')']);
		disp(' ');
		figure;
		subplot(4,1,1);if nDims == 1,plot(z);else imagesc(z);end;xlabel('Data');
		subplot(4,1,2);imagesc(field1);clim([0 1]);xlabel('Raw field');
		subplot(4,1,3);imagesc(field2);clim([0 1]);ylabel(tc);xlabel('Thresholded field');
		subplot(4,1,4);imagesc(field3);clim([0 1]);ylabel(sc);xlabel('Subfield');
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

