function [interpolated,discarded] = Interpolate(samples,timestamps,varargin)

%Interpolate - Interpolate samples (positions, spikes, LFP, etc.) at given timestamps.
%
% Interpolate samples (positions, spikes, LFP, etc.) at given timestamps. As a special
% case, a point process can also be 'interpolated' by changing each original timestamp
% t to the closest target timestamp t* such that t* <= t.
%
%  USAGE
%
%    [interpolated,discarded] = Interpolate(samples,timestamps,<options>)
%
%    samples        a list of samples to interpolate
%    timestamps     a list of timestamps
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'type'        'linear' if samples are linear values (default), or
%                   'circular' otherwise
%     'trim'        discard ('on', default) or keep ('off') samples that would
%                   require extrapolation. If kept, these samples are set to
%                   NaN (this can be useful when the output is required to
%                   have the same size as the input) (default = 'on')
%     'maxGap'      time gaps between original and target times exceeding this
%                   threshold will be ignored (default = Inf)
%    =========================================================================
%
%  OUTPUT
%
%    interpolated   list of interpolated samples
%    discarded      logical vector indicating for which timestamps the samples were
%                   not interpolated (see option 'maxGap')
%
%  SEE
%
%    See also Diff.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
trim = 'on';
maxGap = Inf;
type = 'linear';

if nargin < 2 | mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Interpolate">Interpolate</a>'' for details).');
end
if ~isdvector(timestamps),
	error('Incorrect timestamps - should be a vector (type ''help <a href="matlab:help Interpolate">Interpolate</a>'' for details).');
end
timestamps = timestamps(:);

% Parse parameter list
for j = 1:2:length(varargin),
	if ~ischar(varargin{j}),
		error(['Parameter ' num2str(j+2) ' is not a property (type ''help <a href="matlab:help Interpolate">Interpolate</a>'' for details).']);
	end
	switch(lower(varargin{j})),
		case 'type',
			type = varargin{j+1};
			if ~isstring_FMAT(type,'linear','circular'),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help Interpolate">Interpolate</a>'' for details).');
			end
		case 'trim',
			trim = lower(varargin{j+1});
			if ~isstring_FMAT(trim,'on','off'),
				error('Incorrect value for property ''trim'' (type ''help <a href="matlab:help Interpolate">Interpolate</a>'' for details).');
			end
		case 'maxgap',
			maxGap = varargin{j+1};
			if ~isdscalar(maxGap,'>0'),
				error('Incorrect value for property ''maxGap'' (type ''help <a href="matlab:help Interpolate">Interpolate</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{j}) ''' (type ''help <a href="matlab:help Interpolate">Interpolate</a>'' for details).']);
	end
end

interpolated = [];
if isempty(timestamps) | isempty(samples), return; end

pointProcess = isvector(samples);

% Determine which target timestamps are too far away from sample timestamps (isolated) and should be discarded
[unused,isolated] = Match(timestamps,samples(:,1),'match','closest','error',maxGap);

if pointProcess,
	% 'Interpolate' samples at timestamps (see help above)
	[interpolated,unused,discarded] = Match(samples,timestamps,'error',maxGap);
else
	% Determine which timestamps would require extrapolation
	outside = logical(zeros(size(timestamps)));
	if strcmp(trim,'on'),
		outside = timestamps<samples(1,1)|timestamps>samples(end,1);
	end
	% Discard isolated timestamps and timestamps that would require extrapolation
	discarded = outside|isolated;
	timestamps(discarded) = [];
	% Circular data?
	if strcmp(type,'circular'),
		range = isradians(samples(:,2:end));
		samples(:,2:end) = exp(i*samples(:,2:end));
	end
	% Interpolate samples at timestamps
	interpolated = [timestamps interp1(samples(:,1),samples(:,2:end),timestamps)];
	% Circular data?
	if strcmp(type,'circular'),
		interpolated(:,2:end) = wrap(angle(interpolated(:,2:end)),range);
	end
end