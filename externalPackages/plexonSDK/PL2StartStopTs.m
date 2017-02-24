function [ts] = PL2StartStopTs(filename, startOrStop)
% PL2StartStopTs(filename, startOrStop): read recording start or recording stop event timestamps from a .pl2 file
%
% [ts] = PL2StartStopTs(filename, 'start')
% [ts] = PL2StartStopTs(filename, 'stop')
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   startOrStop - string, either 'start' or 'stop'
%
% OUTPUT:
%   ts - array of timestamps (in seconds)

[ts] = internalPL2StartStopTs(filename, startOrStop);

end