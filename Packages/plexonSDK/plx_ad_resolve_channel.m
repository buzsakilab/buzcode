function [channelNumber] = plx_ad_resolve_chennel(filename, channel)
% plx_ad_resolve_chennel(filename, channel): returns .plx file raw a/d channel number for the specified channel name
%
% [channelNumber] = plx_ad_resolve_chennel(filename, channel)
%
% INPUT:
%   filename - .plx file name
%   channel - 0-based channel number or channel name
%
% OUTPUT:
%   channelNumber - 0-based channel number. 
%                   if channel is a char array (channel name), returns channel number for the specified channel name.
%                      if channel with the specified name is not found, returns -1
%                   if channel is a number, returns this number

channelNumber = -1;

if nargin ~= 2
    error 'expected 2 input arguments';
end

channelNumber = -1;
if ischar(channel) == 1
    [numNames, names] = plx_adchan_names(filename);
    [numMapped, adchans] = plx_ad_chanmap(filename);
    for i=1:numNames
        if strcmp(deblank(names(i,:)), deblank(channel)) == 1
            channelNumber = adchans(i);
            break
        end
    end
else
    channelNumber = channel;
end

end