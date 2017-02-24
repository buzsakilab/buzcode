%Detrend - Detrend signal.
%
%  USAGE
%
%    [detrended,model] = Detrend(samples,<options>)
%
%    samples        list of (t,v1...vn) samples
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'model'       use a predefined model (e.g. from a previous iteration)
%     'window'      size of the regression window (in number of samples)
%     'common'      'on' to use a single model for all windows (default), 'off'
%                   to compute a different model for each window
%     'order'       AR model order (default = 1)
%    =========================================================================
%
%  OUTPUT
%
%    detrended      detrended signal
%    model          AR model (can be reused for other signals)

% Copyright (C) 2012 by Anton Sirota, adapted by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [detrended,model] = Detrend(samples,varargin)

% Defaults
window = [];
model = [];
order = 1;
common = 'on';

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Detrend">Detrend</a>'' for details).');
end

% Check parameter sizes
if ~isdmatrix(samples),
	error('Samples should be a MxN matrix (type ''help <a href="matlab:help Detrend">Detrend</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Detrend">Detrend</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'window',
			window = varargin{i+1};
			if ~isdscalar(window,'>0'),
				error('Incorrect value for property ''window'' (type ''help <a href="matlab:help Detrend">Detrend</a>'' for details).');
			end
		case 'model',
			model = varargin{i+1};
			if ~isdvector(model),
				error('Incorrect value for property ''model'' (type ''help <a href="matlab:help Detrend">Detrend</a>'' for details).');
			end
		case 'order',
			order = varargin{i+1};
			if ~isiscalar(order,'>0'),
				error('Incorrect value for property ''order'' (type ''help <a href="matlab:help Detrend">Detrend</a>'' for details).');
			end
		case 'common',
			common = varargin{i+1};
			if ~isstring_FMAT(common,'on','off'),
				error('Incorrect value for property ''common'' (type ''help <a href="matlab:help Detrend">Detrend</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Detrend">Detrend</a>'' for details).']);
	end
end
%  order = order + 1;

x = samples(:,2:end);
[nSamples,nChannels] = size(x);

y = zeros(nSamples,nChannels);
% List data windows
if isempty(window),
	windows = [1 nSamples];
	nWindows = 1;
else
	nWindows = floor(nSamples/window)+1;
	windows = repmat([1 window],nWindows,1)+repmat([0:nWindows-1]'*window,1,2);
	if nWindows*window > nSamples,
		windows(end,2) = nSamples;
	end   
end

% Process
if ~isempty(model),
	% Model provided as input parameter
	for j = 1:nWindows,
		for i = 1:nChannels,
			y(windows(j,1):windows(j,2),i) = Filter0(model,x(windows(j,1):windows(j,2),i));
		end
	end
elseif strcmp(common,'on'),
	% Use the same model for all windows
	model = arburg(x(windows(1,1):windows(1,2),1),order);
	for j = 1:nWindows,
		for i = 1:nChannels,
			y(windows(j,1):windows(j,2),i) = Filter0(model,x(windows(j,1):windows(j,2),i));
		end
	end
else
	% Determine model for each window
	for j = 1:nWindows,
		for i=1:nChannels,
			model = arburg(x(windows(j,1):windows(j,2),i),order);
			y(windows(j,1):windows(j,2),i) = Filter0(model,x(windows(j,1):windows(j,2),i));
		end
	end
end

detrended = [samples(:,1) y];

% ------------------------------- Helper function -------------------------------

function y = Filter0(b,x)

if mod(length(b),2) ~= 1,
	shift = length(b)/2;
else
	shift = (length(b)-1)/2;
end

[y0,z] = filter(b,1,x);

y = [y0(shift+1:end,:) ; z(1:shift,:)];

