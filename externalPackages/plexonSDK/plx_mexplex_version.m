function [v] = plx_mexplex_version()
% plx_mexplex_version(): returns mexPlex DLL version
%
% [v] = plx_mexplex_version()
%
% INPUT:
%   None
%
% OUTPUT:
%   v - mexPlex DLL version as a number. For example, version 180 corresponds to version 1.8.0.

[v] = mexPlex(29);
