function [AX, AY, A] = AngularAcceleration(VX, VY)


t = Range(VX, 'ts');
ty = Range(VY, 'ts');

if length(ty) ~= length(t)
  error('VX and VY don''t have the same length');
end

if any(ty ~= t)
  error('VX, and VY don''t have the same times.');
end

dt = diff(t);
vx = Data(VX);
vy = Data(VY);

ax = diff(vx) ./ dt;
ax = ax * 10000;
ay = diff(vy) ./ dt;
ay = ay * 10000;

vx = vx(1:end-1);
vy = vy(1:end-1);
v2 = sqrt(vx .* vx + vy .* vy);

ang = (vx .* ay - vy .* ax) ./ (v2 * v2);

A = tsd(t(1:end-1), ang);
AX = tsd(t(1:end-1), ax);
AY = tsd(t(1:end-1), ay);