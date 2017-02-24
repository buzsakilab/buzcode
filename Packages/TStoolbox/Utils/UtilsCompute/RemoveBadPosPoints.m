function [X1, Y1, ds] = RemoveBadPosPoints(X, Y, jump_thresh_x, jump_thresh_y)
% [X1, Y1] = RemoveBadPosPoints(X, Y, jump_thresh)
% removes bad points from Position tsd-s
%
% it removes all the points where a "jump" greater than jump_thresh
% occurs
% INPUTS:
% X, Y: x and y coordinates tsd's (dirty)
% jump_thresh: threshold for jump
% 
% OUTPUTS:
% X1, Y1: x and y coordinates tsd's (clean)

% batta 2001 
% status: under construction
  
  
max_jump = 100;
x = Data(X);
y = Data(Y);
t = Range(X, 'ts');

good_points = RemovePosSkips(x, y, t, jump_thresh_x, jump_thresh_y, 5000);

% $$$ while 1
% $$$   dx = diff(x);
% $$$   dy = diff(y);
% $$$   lx = length(x);
% $$$   dt = diff(t);
% $$$   dx = abs(dx) ./ dt;
% $$$   dy = abs(dy) ./ dt;
% $$$   keyboard
% $$$   
% $$$   good_points = find((dx < jump_thresh_x & dy < jump_thresh_y) | ...
% $$$ 		     dt > 170) + 1;
% $$$   good_points(find(good_points > length(t))) = length(t);
% $$$   x = x(good_points);
% $$$   y = y(good_points);
% $$$   t = t(good_points);
% $$$   
% $$$   keyboard
% $$$   if length(x) >= lx - max_jump
% $$$     break
% $$$   end
% $$$ end

 x = x(good_points);
 y = y(good_points);
 t = t(good_points);

X1 = tsd(t, x);
Y1 = tsd(t, y);