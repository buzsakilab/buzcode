function [errCode] = ddt_write_v(filename, nch, npoints, freq, d)
% ddt_write_v(filename, nch, npoints, freq, d) Write data to a .ddt file
%
% [errCode] = ddt_write_v(filename, nch, npoints, freq, d)
%
% INPUT:
%    filename - if empty string, will use File Open dialog
%    nch - number of channels
%    npoints - number of data points per channel
%    freq - data frequency in Hz
%    d - [nch npoints] data array (in mV)
%
% OUTPUT:
%    errCode - error code: 1 for success, 0 for failure

if nargin < 5
    error 'Expected 5 input arguments';
end
if (isempty(filename))
    [fname, pathname] = uiputfile('*.ddt', 'Specify .ddt file');
    if isequal(fname,0)
        error 'No file was specified'
    end
    filename = fullfile(pathname, fname);
end

[errCode] = mexPlex(25, filename, nch, npoints, freq, d)
