function y = hsv2hsl(x)

%hsv2hsl - Convert hue-saturation-value colors to hue-saturation-luminance.
%
%  USAGE
%
%    y = hsv2hsl(x)
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
  error('Incorrect number of parameters (type ''help <a href="matlab:help hsv2hsl">hsv2hsl</a>'' for details).');
end
if ~isdmatrix(x) || size(x,2) ~= 3 || any(x(:)<0) || any(x(:)>1),
  error('Incorrect HSL matrix (type ''help <a href="matlab:help hsv2hsl">hsv2hsl</a>'' for details).');
end

% Convert HSV to HSL
h = x(:,1);
s = x(:,2);
v = x(:,3);

H = h;
S = s .* v;
L = (2-s).*v;
in = L <= 1;
S(in) = S(in) ./ L(in);
S(in&v==0) = 0;
S(~in) = S(~in) ./ (2-L(~in));
S(~in&s==0&v==1) = 0;
L = L/2;

y = Clip([H S L],0,1);
