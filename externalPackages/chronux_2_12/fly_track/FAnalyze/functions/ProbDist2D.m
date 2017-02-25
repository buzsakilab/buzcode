function [P bins] = ProbDist2D(x, nbins);

%PROBDIST2D
%   Usage: [P bins] = [P bins] = ProbDist2D(x, nbins);
%
% This function calculates the discrete probability distribution for a
% variable that is in a two dimensional phase space.  Because the variable
% exists in a two dimensional phase space, we have to bin in x^2
% (equivalent to dividing P by x). The function simply calculates a
% normalized histogram using 'hist', so nbins can take any form that hist
% can take, although it does not allow for all the other hist options.

%Written by Dan Valente
%September 2007

[H bins] = hist(x.^2, nbins);
N = sum(sum(H));
binsize = bins(3)-bins(2);
P = H./(N*binsize);

%now transform bins back to the unsquared locations.  This is really only
%for conviencence, by no means is it necessary.  If this line is commented
%out, then the plot is simply P vs. x^2 instead of P vs. x
bins = sqrt(bins);  

return;