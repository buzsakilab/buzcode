function [isBuzcode, path] = bz_isBuzcode(varargin)
% USAGE
% [isBuzcode] = bz_isBuzcode(event)
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
if bz_isSession % only worth continuing if it's a session..
    sessionInfo = bz_getSessionInfo;

    if bz_isSessionInfo(sessionInfo) 
        isBuzcode = 1;  % we should also add checks for dat/clu/res/wav files..
    else
        isBuzcode = 0;
    end
else
    isBuzcode = 0;
end