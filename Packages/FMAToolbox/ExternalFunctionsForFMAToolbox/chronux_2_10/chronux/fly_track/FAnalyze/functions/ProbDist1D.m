function [P bins] = ProbDist1D(x, nbins);

%PROBDIST1D
%   Usage: [P bins] = ProbDist1D(x, nbins);
%
% This function calculates the discrete probability distribution for a
% variable that is in a one dimensional phase space.  It is simply a
% normalized histogram using 'hist', so nbins can take any form that hist
% can take, although it does not allow for all the other hist options.

%Written by Dan Valente
%September 2007

[H bins] = hist(x, nbins);
N = sum(sum(H));
binsize = bins(3)-bins(2);
P = H./(N*binsize);  %divide by total points and binsize to normalize correctly

return;