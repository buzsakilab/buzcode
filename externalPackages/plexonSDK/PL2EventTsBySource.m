function [event] = PL2EventTsBySource(filename, source, channelInSource)
% PL2EventTsBySource(filename, source, channelInSource): read event timestamps from a .pl2 file
%                       
% [event] = PL2EventTs(filename, source, channelInSource)
% [event] = PL2EventTs(filename, "KBD', 1)
% [event] = PL2EventTs(filename, 9, 1)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   source - event source name or source numeric id
%   channel -  1-based channel number within source
%           To print the list of event channels, use these commands:
%
%           pl2 = PL2GetFileIndex(filename);
%           PL2Print(pl2.EventChannels);
%
% OUTPUT:
%   event.Ts - array of timestamps (in seconds)
%   event.Strobed - array of strobed event values 

event.Ts = [];
event.Strobed = [];

if nargin ~= 3
    error 'expected 3 input arguments';
end

filename = internalPL2ResolveFilename(filename);
pl2 = PL2GetFileIndex(filename);
channelNumber = internalPL2ResolveChannelBySource(pl2.EventChannels, source, channelInSource);
event = PL2EventTs(filename, channelNumber);

end