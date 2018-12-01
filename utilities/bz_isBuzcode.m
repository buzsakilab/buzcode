function [isBuzcode, path] = bz_isBuzcode(basePath)
% USAGE
% [isBuzcode] = bz_isBuzcode(basePath)
% 
% INPUT
%       basePath    the basePath to check (default: current directory)
%
% OUTPUT
%      
%
% written by david tingley, 2017

if ~exist('basePath','var')
    basePath = pwd;
end

if bz_isSession(basePath) % only worth continuing if it's a session..
    sessionInfo = bz_getSessionInfo(basePath,'noprompts',true);

    if bz_isSessionInfo(sessionInfo) 
        isBuzcode = 1;  % we should also add checks for dat/clu/res/wav files..
    else
        isBuzcode = 0;
    end
else
    isBuzcode = 0;
end
