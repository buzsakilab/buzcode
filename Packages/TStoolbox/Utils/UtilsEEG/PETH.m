function [M, SM, B, center, R] = PETH(X, Y)


binsize = 40 * 10;
nbins = 160 + 1;

M = zeros(nbins, 1);
SM = zeros(nbins, 1);
B = -(nbins - 1) * binsize / 2: binsize:  (nbins - 1) * binsize / 2;
B = B / 10;
R = zeros(nbins, length(Y));

for i = 1: length(Y)
  edges = Y(i)- nbins * binsize / 2: binsize : Y(i) + nbins * binsize / 2;
  Q = histc(X, edges);
  M = M + Q(1:end-1);
  R(:, i) = Q(1:end-1);
  SM = SM + Q(1:end-1) .* Q(1:end-1);
  center(i) = sum(M(((nbins-1)/2)-2:((nbins-1)/2)+2));
end

SM = sqrt(SM - M .* M);

M = M *10000/ (length(Y) * binsize);
SM = SM *10000/ (length(Y) * binsize);
center = center * 10000 / (length(Y) * 5 * binsize);
