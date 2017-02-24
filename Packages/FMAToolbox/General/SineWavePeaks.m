function [t,peaks] = SineWavePeaks(sineWave,varargin)

%SineWavePeaks - Find peaks (or troughs) in a sine wave.
%
% Find the peaks (or troughs) in a sine wave. The algorithm can either determine
% the mid-points between the zero-crossings, or the zeros of the derivative of
% the signal. This assumes a minimum of 10 samples per positive/negative phase.
%
%  USAGE
%
%    [t,peaks] = SineWavePeaks(sineWave,<options>)
%
%    sineWave       an Nx2 matrix of (timestamp,value) pairs
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'mode'        'peaks' or 'troughs' (default = 'peaks')
%     'method'      'zero' for zero-crossings (default), 'diff' for derivative
%     'smooth'      smoothing kernel (in samples) for derivative (default = 5)
%    =========================================================================

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
mode = 'peaks';
method = 'zero';
smooth = 5;

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help SineWavePeaks">SineWavePeaks</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help SineWavePeaks">SineWavePeaks</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'mode',
			mode = lower(varargin{i+1});
			if ~isstring_FMAT(mode,'peaks','troughs'),
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help SineWavePeaks">SineWavePeaks</a>'' for details).');
			end

		case 'method',
			method = lower(varargin{i+1});
			if ~isstring_FMAT(method,'zero','diff'),
				error('Incorrect value for property ''method'' (type ''help <a href="matlab:help SineWavePeaks">SineWavePeaks</a>'' for details).');
			end

		case 'smooth',
			smooth = lower(varargin{i+1});
			if ~isdscalar(smooth),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help SineWavePeaks">SineWavePeaks</a>'' for details).');
			end

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SineWavePeaks">SineWavePeaks</a>'' for details).']);
  end
end


if strcmp(method,'zero'),

	% Find downward and upward going zero-crossings
	previous = sineWave(1:end-1,2);
	current = sineWave(2:end,2);
	down = find(previous > 0 & current < 0);
	up = find(previous < 0 & current > 0);
	% Keep only appropriately ordered pairs of zero-crossings
	if strcmp(mode,'troughs'),
		if up(1) < down(1), up = up(2:end); end
		if down(end) > up(end), down = down(1:end-1); end
	else
		if down(1) < up(1), down = down(2:end); end
		if up(end) > down(end), up = up(1:end-1); end
	end
	t = (sineWave(down,1)+sineWave(up,1))/2;
	peaks = interp1(sineWave(:,1),sineWave(:,2),t);

elseif strcmp(method,'diff'),

	d = diff(Smooth(sineWave(:,2),smooth));
%  	[b,a] = cheby2(4,20,[4 20]*(sineWave(2,1)-sineWave(1,1)));
%  	d = [diff(sineWave(:,1))+sineWave(1:end-1,1) filtfilt(b,a,d)];
	n = length(d);
	D = zeros(n,1);
	for i = 1:10,
%  		disp(['i=' int2str(i) ', D(' int2str(max([i-2 1])) ':' int2str(min([n n+i-3])) ')+=d(' int2str(max([5-i-1 1])) ':' int2str(min([n n-i+3])) ') [lengths: ' int2str(length(max([i-2 1]):min([n n+i-3]))) ', ' int2str(length(max([5-i-1 1]):min([n n-i+3]))) ']']);
		D(max([i-2 1]):min([n n+i-3])) = D(max([i-2 1]):min([n n+i-3])) + d(max([5-i-1 1]):min([n n-i+3]));
	end
	d = [diff(sineWave(:,1))+sineWave(1:end-1,1) D];

	% Find downward and upward going zero-crossings
	previous = d(1:end-1,2);
	current = d(2:end,2);
	if strcmp(mode,'troughs'),
		up = find(previous < 0 & current > 0);
		t = d(up,1);
	else
		down = find(previous > 0 & current < 0);
		t = d(down,1);
	end
	peaks = interp1(sineWave(:,1),sineWave(:,2),t);

%  	figure;hold on;
%  	plot(sineWave(:,1),sineWave(:,2),'b');
%  	plot(d(:,1),d(:,2),'r');
%  	plot(t,zeros(size(t)),'k+','linestyle','none');
%  	for i = 1:length(peaks),
%  		plot([t(i) t(i)],[peaks(i) peaks(i)],'g+','linestyle','none');
%  	end
end