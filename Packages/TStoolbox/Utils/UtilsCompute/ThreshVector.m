function [istart, iend] = ThreshVector(V, threshold, interval, pos)


if pos
  c1 = find(V > threshold);
else
  c1 = find(V < threshold);
end

dc1 = diff(c1);

istart = c1(find(dc1 > interval)+1);
istart = [c1(1); istart];

iend = c1(find(dc1 > interval));
iend = [iend; c1(end)];

ds = iend - istart;
ds = find(ds > interval);
istart = istart(ds);
iend = iend(ds);
