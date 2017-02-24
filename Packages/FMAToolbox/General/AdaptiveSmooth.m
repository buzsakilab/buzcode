function smoothed = AdaptiveSmooth(data,smooth)

%AdaptiveSmooth - Smooth using an adaptive kernel.
%
%  USAGE
%
%    smoothed = AdaptiveSmooth(data,smooth)
%
%    data           data to smooth
%    smooth         vertical and horizontal standard deviations [Sv Sh]
%                   for Gaussian kernel, measured in number of samples
%                   (0 = no smoothing)

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%  vector = min(size(data)) == 1;
vector = isvector(data);
matrix = (~vector & length(size(data)) == 2);
if ~vector & ~matrix,
	error('Smoothing applies only to vectors or matrices (type ''help <a href="matlab:help AdaptiveSmooth">AdaptiveSmooth</a>'' for details).');
end

% Vectors must be 'vertical'
if size(data,1) == 1,
	data = data';
end

if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help AdaptiveSmooth">AdaptiveSmooth</a>'' for details).');
end

if ~isdvector(smooth,'>=0') | (matrix & length(smooth) > 2) | (vector & length(smooth) ~= 1),
	error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help AdaptiveSmooth">AdaptiveSmooth</a>'' for details).');
end

% If Sh = Sv = 0, no smoothing required
if all(smooth==0),
	smoothed = data;
	return
end

if length(smooth) == 1,
	% For 2D data, providing only one value S for the std is interpreted as Sh = Sv = S
	smooth = [smooth smooth];
end

% Build Gaussian kernels
[vSize,hSize] = size(data);
% 1) Vertical kernel
vKernelSize = min([vSize 1001]);
r = (-vKernelSize:vKernelSize)'/vKernelSize;
vKernelStdev = smooth(1)/vKernelSize;
vKernel = exp(-r.^2/(vKernelStdev+eps)^2/2);
vKernel = vKernel/sum(vKernel);
% 2) Horizontal kernel
hKernelSize = min([hSize 1001]);
r = (-hKernelSize:hKernelSize)/hKernelSize;
hKernelStdev = smooth(2)/hKernelSize;
hKernel = exp(-r.^2/(hKernelStdev+eps)^2/2);
hKernel = hKernel/sum(hKernel);
if vector,
	% Vector smoothing
	% Prepend/append data to limit edge effects
	top = flipud(data(1:vKernelSize));
	bottom = flipud(data(end-vKernelSize+1:end));
	data = [top;data;bottom];
	% Convolve (and return central part)
	tmp = conv(vKernel,data);
	n = size(tmp,1);
	d = n - vSize;
	start = d/2+1;
	stop = start + vSize - 1;
	smoothed = tmp(start:stop,:);
else
	% Matrix smoothing
	% Convolve
	if smooth(1) == 0,
		% Smooth only across columns (Sv = 0)
		% Prepend/append data to limit edge effects
		left = fliplr(data(:,1:hKernelSize));
		right = fliplr(data(:,end-hKernelSize+1:end));
		data = [left data right];
		for i = 1:size(data,1),
			tmp = conv(hKernel,data(i,:));
			n = size(tmp,2);
			d = n - hSize;
			start = d/2+1;
			stop = start + hSize - 1;
			smoothed(i,:) = tmp(:,start:stop);
		end
	elseif smooth(2) == 0,
		% Smooth only across lines (Sh = 0)
		% Prepend/append data to limit edge effects
		top = flipud(data(1:vKernelSize,:));
		bottom = flipud(data(end-vKernelSize+1:end,:));
		data = [top;data;bottom];
		for i = 1:size(data,2),
			tmp = conv(vKernel,data(:,i));
			n = size(tmp,1);
			d = n - vSize;
			start = d/2+1;
			stop = start + vSize - 1;
			smoothed(:,i) = tmp(start:stop);
		end
	else
		% Smooth in 2D
		% Prepend/append data to limit edge effects
		top = flipud(data(1:vKernelSize,:));
		bottom = flipud(data(end-vKernelSize+1:end,:));
		data = [top;data;bottom];
		left = fliplr(data(:,1:hKernelSize));
		right = fliplr(data(:,end-hKernelSize+1:end));
		data = [left data right];
		tmp = conv2(vKernel,hKernel,data,'same');
		n = size(tmp,1);
		d = n - vSize;
		vStart = d/2+1;
		vStop = vStart + vSize - 1;
		n = size(tmp,2);
		d = n - hSize;
		hStart = d/2+1;
		hStop = hStart + hSize - 1;
		smoothed = tmp(vStart:vStop,hStart:hStop);
	end
end
