function [isBuzcode] = bz_isBuzcode(events)
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


sessionInfo = bz_getSessionInfo;

if bz_isSessionInfo(sessionInfo) 
    isBuzcode = 1;  % we should also add checks for dat/clu/res/wav files..
else
    isBuzcode = 0;
end