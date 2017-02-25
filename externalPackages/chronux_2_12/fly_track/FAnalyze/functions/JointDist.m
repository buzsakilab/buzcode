function [P bins] = JointDist(x,y,nbins)

%JOINTDIST
%   Usage: [P bins] = JointDist(x,y,nbins)
%
%This function calculates the joint distribution of the variables x and y.
%It is simply a normalized histogram using hist3. nbins can take any form
%of bins that hist3 can take, but no other hist3 options.

%Written by Dan Valente
%October 2007

%ensure x and y are both column matrices
[nrowsx,ncolsx] = size(x);
if (ncolsx ~= 1)
    x = x';
end
[nrowsy,ncolsy] = size(y);
if (ncolsy ~= 1)
    y = y';
end


[H bins] = hist3([x y], nbins);
N = sum(sum(H));
binsize_x = bins{1}(3)-bins{1}(2);
binsize_y = bins{2}(3)-bins{2}(2);
P = H./(N*binsize_x*binsize_y);  %Normalize correctly


return;

