function [event] = PL2EventTs(filename, channel)
% PL2EventTs(filename, channel): read event timestamps from a .pl2 file
%                       
% [event] = PL2EventTs(filename, channel)
% [event] = PL2EventTs(filename, 1)
% [event] = PL2EventTs(filename, 'KBD2')
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - event channel name or 1-based channel number
%           1-based channel number is the index of the event channel in pl2.EventChannels cell array
%           shown in the Name column of the output produced by PL2Print(pl2.EventChannels)
%           To print the list of event channels, use these commands:
%
%           pl2 = PL2GetFileIndex(filename);
%           PL2Print(pl2.EventChannels);
%
% OUTPUT:
%   event.Ts - array of timestamps (in seconds)
%   event.Strobed - array of strobed event values 

event = internalPL2EventTs(filename, channel);

end