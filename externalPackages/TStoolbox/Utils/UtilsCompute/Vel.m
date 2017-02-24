function VX = Vel2d(X)

t = Range(X, 'ts');

dt = diff(t);

vx = diff(Data(X)) ./ dt;

VX = tsd(t(1:end-1), vx);
