function [isSession, path] = bz_isSession(varargin)
% USAGE
% [isSession] = bz_isSession(event)
% 
% checks that the path is a recording session
%
% INPUT
%       
%
% OUTPUT
%      
%
% written by david tingley, 2017

if nargin == 0
    path = pwd;
else
    path = varargin{1};
end
pathSplit = strsplit(path,'/');
session = pathSplit{end}; % get session name 

if exist([lower(session) '.dat']) | exist([lower(session) '.eeg']) | exist([lower(session) '.lfp'])  ...
        | exist([upper(session) '.dat']) | exist([upper(session) '.eeg']) | exist([upper(session) '.lfp']) 
    isSession = 1;  
else
    isSession = 0;
end