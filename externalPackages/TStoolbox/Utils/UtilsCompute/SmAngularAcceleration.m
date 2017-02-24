function [AX, AY, A] = SmAngularAcceleration(VX, VY)


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


V2 = tsd(t(1:end-1), sqrt(vx .* vx + vy .* vy));

V2 = SmoothTsd(V2, 100);
v2 = Data(V2);

VAX = SmoothTsd(tsd(t,vx), 100);
VAY = SmoothTsd(tsd(t,vy), 100);

vx = Data(VX);
vy = Data(VY);
vx = vx(1:end-1);
vy = vy(1:end-1);

% ang = (vx .* ay - vy .* ax) ./ (v2 .* v2);
ang = (vx .* ay - vy .* ax) ./ (v2);
A = tsd(t(1:end-1), ang);
AX = tsd(t(1:end-1), ax);
AY = tsd(t(1:end-1), ay);

A = SmoothTsd(A, 100);
AX = SmoothTsd(AX, 100);
AY = SmoothTsd(AY, 100);
