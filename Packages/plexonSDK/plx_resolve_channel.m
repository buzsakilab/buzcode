function [channelNumber] = plx_resolve_channel(filename, channel)
% plx_resolve_channel(filename, channel): returns .plx file spike channel number for the specified channel name
%
% [channelNumber] = plx_resolve_channel(filename, channel)
%
% INPUT:
%   filename - file name
%   channel - 1-based channel number or channel name
%
% OUTPUT:
%   channelNumber - 1-based channel number. 
%                   if channel is a channel name, returns channel number for the specified channel name.
%                      if channel with the specified name is not found, returns -1
%                   if channel is a number, returns this number

channelNumber = -1;

if nargin ~= 2
    error 'Expected 2 input arguments';
end

channelNumber = -1;
if ischar(channel) == 1
    [numNames, names] = plx_chan_names(filename);
    [numMapped, chans] = plx_chanmap(filename);
    for i=1:numNames
        if strcmp(deblank(names(i,:)), deblank(channel)) == 1
            channelNumber = chans(i);
            break
        end
    end
else
    channelNumber = channel;
end

end

 