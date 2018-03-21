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

if exist([session '.dat']) & exist([session '.eeg']) | exist([session '.lfp'])
    isSession = 1;  
else
    isSession = 0;
end