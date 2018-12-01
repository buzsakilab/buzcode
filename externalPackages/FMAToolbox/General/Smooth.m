function smoothed = Smooth(data,smooth,varargin)

%Smooth - Smooth using a Gaussian kernel.
%
%  USAGE
%
%    smoothed = Smooth(data,smooth,<options>)
%
%    data           data to smooth
%    smooth         vertical and horizontal kernel size:
%                    - gaussian: standard deviations [Sv Sh]
%                    - rect/triangular: window half size [Nv Nh]
%                   (in number of samples, 0 = no smoothing)
%    <options>      optional list of property-value pairs (see table below))
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'type'        two letters (one for X and one for Y) indicating which
%                   coordinates are linear ('l') and which are circular ('c')
%                   - for 1D data, only one letter is used (default 'll')
%  	'kernel'		  either 'gaussian' (default), 'rect' (running average), 
%  					  or 'triangular' (weighted running average)
%    =========================================================================
%

% Copyright (C) 2004-2015 by Michaël Zugaro, 2013 Nicolas Maingret
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

maxSize = 10001;

if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
end

vector = isvector(data);
matrix = (~vector & length(size(data)) == 2);
if ~vector & ~matrix,
	error('Smoothing applies only to vectors or matrices (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
end

% Vectors must be 'vertical'
if size(data,1) == 1,
	data = data';
end

% Default values
kernel = 'gaussian';
if vector, type = 'l'; else type = 'll'; end

if ~isdvector(smooth,'>=0') | (matrix & length(smooth) > 2) | (vector & length(smooth) ~= 1),
	error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
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

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'type',
			type = lower(varargin{i+1});
			if (vector && ~isstring_FMAT(type,'c','l')) || (~vector && ~isstring_FMAT(type,'cc','cl','lc','ll')),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
			end
		case 'kernel',
			kernel = lower(varargin{i+1});
			if ~isstring_FMAT(kernel,'gaussian','rectangular','triangular'),
				error('Incorrect value for property ''kernel'' (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Smooth">Smooth</a>'' for details).']);

  end
end

% Build kernels
[vSize,hSize] = size(data);

if strcmp(kernel,'rectangular'),

	% Rectangular kernel (running average)
	% 1) Vertical kernel
	vKernelSize = 2*smooth(1)+1;
	vKernel = ones(vKernelSize,1);
	vKernel = vKernel/sum(vKernel);
	% 2) Horizontal kernel
	if ~vector,
		hKernelSize = 2*smooth(2)+1;
		hKernel = ones(hKernelSize,1);
		hKernel = hKernel/sum(hKernel);
	end

elseif strcmp(kernel,'triangular'),

	% Triangular kernel
	% 1) Vertical kernel
	vKernelSize = 2*smooth(1)+1;
	vKernel = triang(vKernelSize);
	vKernel = vKernel/sum(vKernel);
	% 2) Horizontal kernel
	if ~vector,
		hKernelSize = 2*smooth(2)+1;
		hKernel = ones(hKernelSize,1);
		hKernel = hKernel/sum(hKernel);
	end

else

	% Gaussian kernel
	% 1) Vertical kernel
	if vSize > maxSize, warning(['Kernel too large; using ' int2str(maxSize) ' points.']); end
	vKernelSize = min([vSize maxSize]);
	r = (-vKernelSize:vKernelSize)'/vKernelSize;
	vKernelStdev = smooth(1)/vKernelSize;
	vKernel = exp(-r.^2/(vKernelStdev+eps)^2/2);
	vKernel = vKernel/sum(vKernel);
	% 2) Horizontal kernel
	if ~vector,
		if hSize > maxSize, warning(['Kernel too large; using ' int2str(maxSize) ' points.']); end
		hKernelSize = min([hSize maxSize]);
		r = (-hKernelSize:hKernelSize)/hKernelSize;
		hKernelStdev = smooth(2)/hKernelSize;
		hKernel = exp(-r.^2/(hKernelStdev+eps)^2/2);
		hKernel = hKernel/sum(hKernel);
	end
	
end

if vector,
	% Vector smoothing
	% Prepend/append data to limit edge effects
    %% This changes the behavior of Smooth, to avoid linear variable edge effects
	if strcmp(type,'l'),
		% For linear data, flip edge data
		% top = 2*data(1)-flipud(data(1:vKernelSize));
		% bottom = 2*data(end)-flipud(data(end-vKernelSize+1:end));
% 		top = flipud(data(1:vKernelSize));
% 		bottom = flipud(data(end-vKernelSize+1:end));
	else
		% For circular data, wrap edge data
		top = data(end-vKernelSize+1:end);
		bottom = data(1:vKernelSize);
        data = [top;data;bottom];
	end
% 	data = [top;data;bottom];
%%
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
		if strcmp(type(1),'l'),
			% For linear data, flip edge data
			% left = 2*repmat(data(:,1),1,hKernelSize)-fliplr(data(:,1:hKernelSize));
			% right = 2*repmat(data(:,end),1,hKernelSize)-fliplr(data(:,end-hKernelSize+1:end));
			left = fliplr(data(:,1:hKernelSize));
			right = fliplr(data(:,end-hKernelSize+1:end));
		else
			% For circular data, wrap edge data
			left = data(:,end-hKernelSize+1:end);
			right = data(:,1:hKernelSize);
		end
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
		if strcmp(type(2),'l'),
			% For linear data, flip edge data
			% top = 2*repmat(2*data(1,:),vKernelSize,1)-flipud(data(1:vKernelSize,:));
			% bottom = 2*repmat(2*data(end,:),vKernelSize,1)-flipud(data(end-vKernelSize+1:end,:));
			top = flipud(data(1:vKernelSize,:));
			bottom = flipud(data(end-vKernelSize+1:end,:));
		else
			% For circular data, wrap edge data
			bottom = data(1:vKernelSize,:);
			top = data(end-vKernelSize+1:end,:);
		end
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
		if strcmp(type(2),'l'),
			% For linear data, flip edge data
  			% top = 2*repmat(2*data(1,:),vKernelSize,1)-flipud(data(1:vKernelSize,:));
  			% bottom = 2*repmat(2*data(end,:),vKernelSize,1)-flipud(data(end-vKernelSize+1:end,:));
			top = flipud(data(1:vKernelSize,:));
			bottom = flipud(data(end-vKernelSize+1:end,:));
		else
			% For circular data, wrap edge data
			bottom = data(1:vKernelSize,:);
			top = data(end-vKernelSize+1:end,:);
		end
		data = [top;data;bottom];
		if strcmp(type(1),'l'),
			% For linear data, flip edge data
			% left = 2*repmat(2*data(:,1),1,hKernelSize)-fliplr(data(:,1:hKernelSize));
  			% right = 2*repmat(2*data(:,end),1,hKernelSize)-fliplr(data(:,end-hKernelSize+1:end));
			left = fliplr(data(:,1:hKernelSize));
			right = fliplr(data(:,end-hKernelSize+1:end));
		else
			% For circular data, wrap edge data
			left = data(:,end-hKernelSize+1:end);
			right = data(:,1:hKernelSize);
		end
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
