function [VX, VY] = Vel2d(X, Y)

t = Range(X, 'ts');
ty = Range(Y, 'ts');

if length(ty) ~= length(t)
  error('X and Y don''t have the same length');
end

if any(ty ~= t)
  error('X, and Y don''t have the same times.');
end

dt = diff(t);

vx = diff(Data(X)) ./ dt;
vy = diff(Data(Y)) ./ dt;

VX = tsd(t(1:end-1), vx);
VY = tsd(t(1:end-1), vy);
