function [nch, npoints, freq, d] = ddt(filename)
% ddt(filename) Read data from a .ddt file
%
% [nch, npoints, freq, d] = ddt(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   nch - number of channels
%   npoints - number of data points for each channel
%   freq - A/D frequency
%   d - [nch npoints] data array 

if nargin < 1
    error 'Expected 1 input argument';
end
if (isempty(filename))
   [fname, pathname] = uigetfile('*.ddt', 'Select a Plexon .ddt file');
   if isequal(fname,0)
     error 'No file was selected'
   end
   filename = fullfile(pathname, fname);
end

[nch, npoints, freq, d] = mexPlex(1, filename);