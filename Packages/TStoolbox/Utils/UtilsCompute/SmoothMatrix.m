function S = SmoothMatrix(V, l)

l = floor(l / 2) * 2;

hh = hamming(l);

mv = Data(V);

t = Range(V, 'ts');
mvs = zeros(size(mv));

for i = 1:size(mv,2)
  v = mv(:,i);
  v = conv(v, hh);

  v = v(l/2:end-l/2);

  v = v / sum(hh);

  mvs(:,i) = v;
end

S = tsd(t, mvs);
