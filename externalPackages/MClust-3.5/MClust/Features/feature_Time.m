function [TimeData, TimeNames,TimePars] = feature_Time(V, ttChannelValidity, Params)

% MClust
% [TimeData, TimeNames] = feature_Peak(V, ttChannelValidity)
% Calculate time feature 
%
% INPUTS
%    V = TT tsd
%    ttChannelValidity = nCh x 1 of booleans
%
% OUTPUTS
%    Data - nSpikes 
%    Names - "Time"
%
% ADR April 2008
% version M1.0
% RELEASED as part of MClust 3.5
% See standard disclaimer in Contents.m

TimeData = Range(V, 'sec');
TimeNames = {'_Time'};
TimePars = [];