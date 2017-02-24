function [n, ts] = plx_ts(filename, channel, unit)
% plx_ts(filename, channel, unit): read spike timestamps from a .plx file
%
% [n, ts] = plx_ts(filename, channel, unit)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - 1-based channel number or channel name
%   unit  - unit number (0- unsorted, 1-4 units a-d)
%
% OUTPUT:
%   n - number of timestamps
%   ts - array of timestamps (in seconds)

n = 0;
ts = -1;

if nargin ~= 3
    error 'expected 3 input arguments';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    ts = PL2Ts( filename, channel, unit );
    if numel(ts) == 0
        ts = -1;
    else
        n = numel(ts);
    end
    return;
end

channelNumber = plx_resolve_channel(filename, channel);
if channelNumber == -1
    fprintf('\n plx_ts: no header for the specified spike channel.\n');
    return
end

[n, ts] = mexPlex(5, filename, channelNumber, unit);

end