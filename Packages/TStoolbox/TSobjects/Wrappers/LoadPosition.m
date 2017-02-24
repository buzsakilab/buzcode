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

Fs = 1250/32;

whlt = dlmread([fbasename '.whl']);
% keyboard
[whl GoodRanges] = CleanWhlForR(whlt);

t = (1:size(whlt,1))'/Fs;
GoodRanges = GoodRanges/Fs;
