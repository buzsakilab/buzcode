function y = hsl2hsv(x)

%hsl2hsv - Convert hue-saturation-luminance colors to hue-saturation-value.
%
%  USAGE
%
%    y = hsl2hsv(x)
%
%    x              Nx3 RGB matrix (values in [0..1])

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help hsl2hsv">hsl2hsv</a>'' for details).');
end
if ~isdmatrix(x) || size(x,2) ~= 3 || any(x(:)<0) || any(x(:)>1),
  error('Incorrect HSL matrix (type ''help <a href="matlab:help hsl2hsv">hsl2hsv</a>'' for details).');
end

% Convert from HSL to HSV

h = x(:,1);
s = x(:,2);
l = x(:,3);

H = h;
S = s;
in = l <= 0.5;
S(in) = S(in) .* l(in);
S(~in) = S(~in) .* (1-l(~in));
V = l + S;
in = V > 0;
S(in) = 2 * S(in) ./ V(in);
S(~in) = 0;


y = Clip([H S V],0,1);


%  in = l <= 0.5;
%  s(in) = s(in) .* l(in);
%  s(~in) = s(~in) .* (1-l(~in));
%  l = l + s;
%  
%  in = l > 0;
%  s(in) = 2 * s(in) ./ l(in);
%  s(~in) = 0;
%  
%  H = h;
%  S = s;
%  V = l;
