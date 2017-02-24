function [n] = plx_close(filename)
% plx_close(filename): close the .plx file
%
% [n] = plx_close(filename)
%
% INPUT:
%   filename - if empty string, will close any open files
%
% OUTPUT:
%   n - always 0

[n] = mexPlex(22, filename);
