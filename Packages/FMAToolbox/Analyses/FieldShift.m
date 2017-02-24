function [relativeShift,absoluteShift,xc] = FieldShift(control,test,varargin)

%FieldShift - Estimate firing field shifts between two conditions.
%
%  For each place cell, the absolute field shift between the two conditions is
%  estimated as the mode of the spatial cross-correlogram between the
%  respective firing fields.
%
%  This value is then divided by the average field size across conditions,
%  yielding a relative shift. This measures by how much the field moves between
%  the two conditions, in proportion to the field size.
%
%  Current implementation is only for 1D environments.
%
%  USAGE
%
%    [relativeShift,absoluteShift,xc] = FieldShift(control,test,<options>)
%
%    control        firing fields in control condition (MxN: M fields, N bins)
%    test           firing fields in test condition
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'type'        'linear' for linear tracks (default), 'circular' otherwise
%     'show'        plot results (default = 'off')
%    =========================================================================
%

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
type = 'linear';
show = 'off';

% Check number of parameters
if nargin < 3 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help FieldShift">FieldShift</a>'' for details).');
end

% Check parameter sizes
if ~isdmatrix(control) || ~isdmatrix(test),
	error('Both sets of firing fields should be MxN matrices (type ''help <a href="matlab:help FieldShift">FieldShift</a>'' for details).');
end
if ~all(size(control)==size(test)),
	error('Both sets of firing fields should have the same sizes (type ''help <a href="matlab:help FieldShift">FieldShift</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FieldShift">FieldShift</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'type',
			type = varargin{i+1};
			if ~isstring_FMAT(type,'linear','circular'),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help FieldShift">FieldShift</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring_FMAT(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help FieldShift">FieldShift</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FieldShift">FieldShift</a>'' for details).']);
	end
end

% Compute cross-correlations between fields
[m,n] = size(control);
[xc,mode] = XCorr1(test,control,'type',type); % test vs control yields xc>0 for rightward shift

% Compute mean field sizes across control and repeat, and across control and test
fieldSize = mean([sum(control>0,2) sum(test>0,2)],2)/n;
% Compute field shift as the distance relative to size
relativeShift = mode ./ fieldSize;
absoluteShift = mode;

% Find peak positions (optimized code, uneasy to read)
[y,x] = find(control==repmat(max(control,[],2),1,n));
peakControl = Accumulate(y,x)./Accumulate(y,1)/n;
[y,x] = find(test==repmat(max(test,[],2),1,n));
peakTest = Accumulate(y,x)./Accumulate(y,1)/n;

if strcmp(show,'on'),

	% Plot cross-correlograms
	figure;subplot(1,2,1);hold on;
	PlotColorCurves(xc,[-1 1],0,'w--',mode,'.k');

	% Plot field modes in control vs test (density plot)
	subplot(1,2,2);hold on;
	nBins = 100;
	range = [0 1];
	c = Bin(peakControl,range,nBins);
	t = Bin(peakTest,range,nBins);
	z = Accumulate([c t],1,[nBins nBins]);
	z = Smooth([z z;z z],3)'/n;
	PlotColorMap(z,'x',linspace(range(1),2*range(2),2*n),'y',linspace(range(1),2*range(2),2*n));
	xlabel('Control Peak Position');
	ylabel('Test Peak Position');
	clim([0 0.3*max(z(:))]);
	[beta,R2] = CircularRegression(peakControl*2*pi,peakTest*2*pi,'slope',1);
	hold on;
	PlotSlope(1,1,beta(1),2,'w--');
	
end
