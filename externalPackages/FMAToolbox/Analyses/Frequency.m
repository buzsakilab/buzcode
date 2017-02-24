function frequency = Frequency(timestamps,varargin)

%Frequency - Compute instantaneous frequency for a point process (e.g. a spike train).
%
% Compute instantaneous frequency either by smoothing the point process using a fixed
% or adaptive Gaussian kernel, or by simply taking the inverse of the interevent intervals.
% Adaptive kernel smoothing is described in Richmond et al. (1990).
%
%  USAGE
%
%    frequency = Frequency(timestamps,<options>)
%
%    timestamps     list of timestamps
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'method'      convolute with a fixed ('fixed', default) or adaptive
%                   ('adaptive') Gaussian kernel, or compute inverse inter-
%                   event intervals ('inverse')
%     'limits'      [start stop] in seconds (default = approx. first and last
%                   timestamps)
%     'binSize'     bin size in seconds (default = 0.050)
%     'smooth'      Gaussian kernel width in number of samples (default = 2)
%     'show'        plot results (default = 'off')
%    =========================================================================
%
%  NOTE
%
%    While the convolution method returns a Mx2 matrix (time bins in column 1,
%    frequencies in column 2), the inverse interspike intervals method returns
%    a Nx2 matrix (original spike timestamps in column 1, frequencies in column 2)

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
method = 'fixed';
binSize = 0.05;
smooth = 2;
limits = [];
show = 'off';

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Frequency">Frequency</a>'' for details).');
end

% Check parameter sizes
if size(timestamps,2) ~= 1,
	error('Parameter ''timestamps'' is not a vector (type ''help <a href="matlab:help Frequency">Frequency</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Frequency">Frequency</a>'' for details).']);
	end
	switch(lower(varargin{i})),

		case 'method',
			method = lower(varargin{i+1});
			if ~isstring_FMAT(method,'fixed','adaptive','inverse','iisi'),
				error('Incorrect value for property ''method'' (type ''help <a href="matlab:help Frequency">Frequency</a>'' for details).');
			end

		case 'binsize',
			binSize = varargin{i+1};
			if ~isdscalar(binSize,'>=0'),
				error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help Frequency">Frequency</a>'' for details).');
			end

		case 'limits',
			limits = varargin{i+1};
			if ~isdvector(limits,'#2'),
				error('Incorrect value for property ''limits'' (type ''help <a href="matlab:help Frequency">Frequency</a>'' for details).');
			end

		case 'smooth',
			smooth = varargin{i+1};
			if ~isdscalar(smooth,'>=0'),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help Frequency">Frequency</a>'' for details).');
			end

		case 'show',
			show = varargin{i+1};
			if ~isstring_FMAT(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help Frequency">Frequency</a>'' for details).');
			end

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Frequency">Frequency</a>'' for details).']);

  end
end

if isempty(limits),
	limits = [timestamps(1)-10*binSize timestamps(end)+10*binSize];
end

if strcmp(method,'inverse') | strcmp(method,'iisi'),
	% Inter-spike intervals
	ds = diff(timestamps);
	f1 = 1./ds;
	t = timestamps(1:end-1) + ds/2;
	f2 = Interpolate([t f1],timestamps(2:end-1));
	f2 = [timestamps(1) 0;f2;timestamps(end) 0];
	frequency = f2;
else
	% Smoothing
	% 1) Fixed-kernel
	t = (limits(1):binSize:limits(2))';
	T = Bin(timestamps,t);
	binned = Accumulate(T,1,size(t));
	f = Smooth(binned/binSize,smooth);
	frequency = [t f];
	if strcmp(method,'adaptive'),
		% 2) Variable-kernel (requires the above 'pilot' fixed-kernel estimate)
		% Compute variable-kernel sigma
		N = length(t);
		mu = exp(1/N*sum(log(f(f~=0)))); % Correction: geometric mean is computed using only non-zero values
		lambda = sqrt(f/mu);
		sigma = (smooth*binSize)./lambda;
		sigma = Clip(sigma,0,N*binSize/3);
		% Perform variable-kernel smoothing
		binned = [flipud(binned);binned;flipud(binned)];
		for i = 1:N,
			% Gaussian
			x = (0:binSize:3*sigma(i))';
			x = [flipud(-x);x(2:end)];
			kernel = exp(-x.^2/sigma(i)^2);
			kernel = kernel/sum(kernel);
			% 'Pointwise' convolution
			n = (length(kernel)-1)/2;
			bins = N + i + (-n:n);
			S2(i) = sum(binned(bins).*kernel)/binSize;
		end
		frequency = [t S2'];
	end
end

if strcmp(lower(show),'on'),
	PlotXY(frequency);
	hold on;
	PlotTicks(timestamps,'size',10,'k');
end
