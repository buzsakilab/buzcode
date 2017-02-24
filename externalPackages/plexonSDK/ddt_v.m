function [nch, npoints, freq, d] = ddt_v(filename)
% ddt_v(filename) Read data from a .ddt file returning samples in mV
%
% [nch, npoints, freq, d] = ddt_v(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   nch - number of channels
%   npoints - number of data points for each channel
%   freq - A/D frequency
%   d - [nch npoints] data array (in mV)

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

[nch, npoints, freq, d] = mexPlex(20, filename);