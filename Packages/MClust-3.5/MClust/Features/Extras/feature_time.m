function [timeData, timeNames, timePars] = feature_time(V, ttChannelValidity)

% MClust% [timeData, timeNames] = feature_time(V, ttChannelValidity)% Calculate time feature %% INPUTS%    V = TT tsd%    ttChannelValidity = nCh x 1 of booleans%% OUTPUTS%    Data - nSpikes times %    Names - "Time"
%
% ADR April 1998%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.


timeData = Range(V, 'ts');
timeNames = {'time (ts)'};
timePars = {};
