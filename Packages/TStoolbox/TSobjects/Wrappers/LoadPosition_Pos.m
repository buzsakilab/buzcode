function [whl,t,GoodRanges] = LoadPosition(fbasename)

% Intermediate wrapper, should not be used directly
%
% USAGE
% [whl,t,GoodRanges,ep] = LoadPosition(fbasename)
% 
% OUTPUT:
%   whl: the position matrix
%   t: time vector
%   Good Ranges: 

whlt = dlmread([fbasename '.pos']);
whlt(isnan(whlt)) = -1;
[whl,GoodRanges] = CleanWhlForR(whlt(:,2:end));

t = whlt(:,1);
GoodRanges = [t(GoodRanges(:,1)) t(GoodRanges(:,2))];
