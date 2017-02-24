function m = AdjustColorMap(varargin)

%AdjustColorMap - Adjust colormap for current figure, i.e. change gamma.
%
%  USAGE
%
%    m = AdjustColorMap(<options>)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'shift'       proportional hue shift (0 = no change, default)
%     'gamma'       standard gamma for luminance (1 = no change, default)
%     'hgamma'      gamma for hue (1 = no change, default)
%     'hrotate'     proportional hue rotation (0 = no change, default)
%    =========================================================================

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
gamma = 1;
hgamma = 1;
hrotate = 1;
shift = 0;

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help AdjustColorMap">AdjustColorMap</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'gamma',
			gamma = varargin{i+1};
			if ~isdscalar(gamma,'>=0'),
			error('Incorrect value for property ''gamma'' (type ''help <a href="matlab:help AdjustColorMap">AdjustColorMap</a>'' for details).');
			end
		case 'hgamma',
			hgamma = varargin{i+1};
			if ~isdscalar(hgamma,'>=0'),
			error('Incorrect value for property ''hgamma'' (type ''help <a href="matlab:help AdjustColorMap">AdjustColorMap</a>'' for details).');
			end
		case 'shift',
			shift = varargin{i+1};
			if ~isdscalar(shift),
			error('Incorrect value for property ''shift'' (type ''help <a href="matlab:help AdjustColorMap">AdjustColorMap</a>'' for details).');
			end
		case 'hrotate',
			hrotate = varargin{i+1};
			if ~isdscalar(hrotate),
			error('Incorrect value for property ''hrotate'' (type ''help <a href="matlab:help AdjustColorMap">AdjustColorMap</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help AdjustColorMap">AdjustColorMap</a>'' for details).']);
	end
end

% Convert from RGB to HSL
hsl = hsv2hsl(rgb2hsv(colormap));
h = hsl(:,1);
s = hsl(:,2);
l = hsl(:,3);

% Correct luminance and hue gammas
h = h .^ (1/hgamma);
l = l .^ (1/gamma);

% Convert from HSL to RGB
m = hsv2rgb(hsl2hsv([h s l]));
%  bad = any(isnan(m),2);
%  m(bad,:) = [];
%  m = Clip(m,0,1);

% Rotate colors
n = size(m,1);
r = mod(round(abs(hrotate)*n),n);
m = [m(r+1:end,:) ; m(1:r,:)];

% Correct hue
l = size(m,1);
n = round(abs(shift)*l);
if shift > 0,
	m = [m ; repmat(m(end,:),n,1)];
elseif shift < 0,
	m = [repmat(m(1,:),n,1) ; m];
end

% Apply
colormap(m);








%  % Convert RGB to HSV
%  hsv = rgb2hsv(colormap);
%  
%  % Convert HSV to HSL
%  h = hsv(:,1);
%  s = hsv(:,2) .* hsv(:,3);
%  l = (2-hsv(:,2)).*hsv(:,3);
%  if l <= 1,
%  	s = s ./ l;
%  else
%  	s = s ./ (2-l);
%  end
%  l = l/2;
%  
%  % Correct gamma
%  l = l .^ (1/gamma);
%  
%  % Convert HSL to HSV
%  hsv(:,1) = h .^ (1/hgamma);
%  l = l * 2;
%  if l <= 1,
%  	s = s .* l;
%  else
%  	s = s .* (2-l);
%  end
%  hsv(:,2) = (2 * s) ./ (l + s);
%  hsv(:,3) = (l + s) / 2;
%  
%  % Convert to RGB
%  m = hsv2rgb(hsv);
%  bad = any(isnan(m),2);
%  m(bad,:) = [];
%  m = Clip(m,0,1);
%  
%  % Rotate colors
%  n = size(m,1);
%  r = mod(round(abs(hrotate)*n),n);
%  m = [m(r+1:end,:) ; m(1:r,:)];
%  
%  % Correct hue
%  l = size(m,1);
%  n = round(abs(shift)*l);
%  if shift > 0,
%  	m = [m ; repmat(m(end,:),n,1)];
%  elseif shift < 0,
%  	m = [repmat(m(1,:),n,1) ; m];
%  end
%  
%  % Apply
%  colormap(m);
