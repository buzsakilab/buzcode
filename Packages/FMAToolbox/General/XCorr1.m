function [xc,mode] = XCorr1(x,y,varargin)

%XCorr1 - Compute cross-correlograms and their modes for pairs of curves.
%
%  USAGE
%
%    [xc,mode] = XCorr1(x,y,<options>)
%
%    x              first set of curves (MxN matrix: M curves, N bins)
%    y              second set of curves
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'type'        'linear' for linear curves (default), 'circular' otherwise
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

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help XCorr1">XCorr1</a>'' for details).');
end

% Check parameter sizes
if ~isdmatrix(x) || ~isdmatrix(y),
	error('Both sets of curves should be MxN matrices (type ''help <a href="matlab:help XCorr1">XCorr1</a>'' for details).');
end
if ~all(size(x)==size(y)),
	error('Both sets of curves should have the same sizes (type ''help <a href="matlab:help XCorr1">XCorr1</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help XCorr1">XCorr1</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'type',
			type = varargin{i+1};
			if ~isstring_FMAT(type,'linear','circular'),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help XCorr1">XCorr1</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help XCorr1">XCorr1</a>'' for details).']);
	end
end

[m,n] = size(x);

if strcmp(type,'linear'),
	xc = nan(m,2*n-1);
	for i = 1:m,
		% Cross-correlagram i
		xci = xcorr(x(i,:),y(i,:),'coeff');
		% mode of cross-correlogram i (in [-1..1])
		mi = round(mean(find(xci==max(xci))));
		mi = mi/n-1;
		% Store in matrices
		xc(i,:) = xci;
		mode(i,1) = mi;
	end
else
	xc = nan(m,2*n-1);
	for i = 1:m,
		% Cross-correlagram i
		xci = xcorr(x(i,:),y(i,:),'coeff');
		% mode of cross-correlogram i (in [-1..1])
		mi = round(mean(find(xci==max(xci))));
		mi = mi/n-1;
		% For circular curves, very large modes must be corrected
		if abs(mi) > 0.5,
			if mi > 0,
				xci = [xci(n+1:end) xci(1:n)];
				mi = mi - 1;
			else
				xci = [xci(n:end) xci(1:n-1)];
				mi = mi +1;
			end
		end
		xc(i,:) = xci;
		mode(i,1) = mi;
	end
end


%  	xc = nan(m,2*n-1);
%  	for i = 1:m,
%  		% Cross-correlagram i
%  		X = repmat(x(i,:),n,1);
%  		Y = repmat(y(i,:),1,n+1);
%  		Y(n+1:n+1:end) = [];
%  		Y = reshape(Y,n,[]);
%  		xci = sum(X.*Y,2)';
%  		xci = [xci(2:end) xci];
%  		% mode of cross-correlogram i (in [0..1])
%  		% (we need to compute mean values of x, but x is circular, so we transform it to angles)
%  		a = find(xci==max(xci));
%  		a = exp(j*a*2*pi/n);
%  		mi = round(wrap(angle(mean(a)),2)/(2*pi)*n)/n-1;
%  		if abs(mi) > 0.5,
%  			if mi > 0,
%  				mi = mi - 1;
%  			else
%  				mi = mi +1;
%  			end
%  		end
%  		% Zero outside primary field
%  		m = round((mi+1)*n);
%  		h = ceil(n/2);
%  		xci(m+h:end) = 0;
%  		xci(1:m-h) = 0;
%  		% Store in matrices
%  		xc(i,:) = xci;
%  		mode(i,1) = mi;
%  	end
%  	% Scale so that autocorrelation at shift 0 is 1
%  	x0 = sum(x.^2,2);
%  	y0 = sum(y.^2,2);
%  	scale = sqrt(x0.*y0);
%  	xc = xc./repmat(scale,1,2*n-1);
%  
%  	end
