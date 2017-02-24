function [S, E] = VelFilter(V, threshold)
% [S, E] = VelFilter(V, threshold)
%
% Start and end times of hte epochs in which te speed il less than threshold
%
% INPUTS:
% V = tsd absolute velocity
% threshold = ...
% 
% OUTPUTS:
% S = Column vector of start times
% E = Column vector of end times

% batta 2000 
% status alpha


t = Range(V, 'ts');
v = Data(V);

v1 = [0 v'];
v2 = [v' 0];

s = find(v2 < threshold & v1 > threshold);
e = find(v2 > threshold & v1 < threshold);

if (s(1) == 1)
  s = s(2:length(s));
end

j = 1;
i = 1;
ii = 1;

len_s = length(s)

while i < len_s
  while(e(j) < s(i))  & j < len_s
    j = j+ 1;
  end
  
  if e(j)  < s(i+1)
    S(ii) = s(i);
    E(ii) = e(j);
    ii = ii + 1;
  end
  i = i + 1;
end

if length(E) ~= length(S)
  if length(E) == length(S) - 1
    E(ii-1) = len_s;
  else
    error('S and E lengths differs');
  end
end

S = t(S);
E = t(E);
