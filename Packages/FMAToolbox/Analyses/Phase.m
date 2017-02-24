function [phase,amplitude,unwrapped] = Phase(samples,times)

%Phase - Compute instantaneous phase in signal.
%
% Compute instantaneous phase (wrapped and unwrapped) and amplitude
% in signal.
%
%  USAGE
%
%    [phase,amplitude,unwrapped] = Phase(samples,times)
%
%    samples        signal, e.g. filtered local field potential samples
%    times          optional timestamps where phase should be interpolated
%
%  NOTE
%
%    Angles are returned in radians.
%
%  SEE
%
%    See also FilterLFP, PhasePrecession, PhaseMap.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Phase">Phase</a>'' for details).');
end

% Check parameter sizes
if size(samples,2) < 2,
	error('Parameter ''samples'' is not a matrix (type ''help <a href="matlab:help Phase">Phase</a>'' for details).');
end
if nargin == 2,
	if ~isdvector(times),
		error('Parameter ''times'' is not a vector (type ''help <a href="matlab:help Phase">Phase</a>'' for details).');
	end
	if size(times,2) ~= 1, times = times'; end
end

phase = zeros(size(samples));
phase(:,1) = samples(:,1);
unwrapped = zeros(size(samples));
unwrapped(:,1) = samples(:,1);
amplitude = zeros(size(samples));
amplitude(:,1) = samples(:,1);

for j = 2:size(samples,2),
	% Compute phase and amplitude using Hilbert transform
	h = hilbert(samples(:,j));
	phase(:,j) = mod(angle(h),2*pi);
	amplitude(:,j) = abs(h);
	unwrapped(:,j) = unwrap(phase(:,j));
end

if nargin == 2,
%  	phase = Interpolate([phase(:,1) exp(i*phase(:,2))],times(:,1));
%  	phase(:,2) = angle(phase(:,2));
	phase = Interpolate(phase,times(:,1),'type','circular');
	amplitude = Interpolate(amplitude,times(:,1));
	unwrapped = Interpolate(unwrapped,times(:,1));
end