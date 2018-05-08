function [isSessionInfo] = bz_isSessionInfo(SessionInfo)
% USAGE
% [isSessionInfo] = bz_isSessionInfo(SessionInfo)
% 
% INPUT
%       SessionInfo   - struct with the following fields
%                   .channels
%                   .region
%                   .depth
%                   .nChannels
%                   .FileName
%                   .Date
%                    
%
% OUTPUT
%      logical true if struct meets event criteria, false if otherwise
%
% written by david tingley, 2017


if isfield(SessionInfo,'channels')
     if isvector(SessionInfo.channels) 
         isSessionInfo = true;
     else
         warning('.channels field is not formatted correctly')
         isSessionInfo = false;
     end
end
if isfield(SessionInfo,'region')
     if iscell(SessionInfo.region) 
         isSessionInfo = true;
     else
         warning('.region field is not formatted correctly')
         isSessionInfo = false;
     end
end
if isfield(SessionInfo,'depth')
     if isnumeric(SessionInfo.depth) 
         isSessionInfo = true;
     else
         warning('.depth field is not formatted correctly')
         isSessionInfo = false;
     end
end
if isfield(SessionInfo,'spikeGroups')
     if isstruct(SessionInfo.spikeGroups) % should we go deeper here to check the struct?
         isSessionInfo = true;
     else
         warning('.spikeGroups field is not formatted correctly')
         isSessionInfo = false;
     end
end
if isfield(SessionInfo,'FileName')
     if ischar(SessionInfo.FileName) 
         isSessionInfo = true;
     else
         warning('.FileName field is not formatted correctly')
         isSessionInfo = false;
     end
end
if isfield(SessionInfo,'Date')
     if ischar(SessionInfo.Date) 
         isSessionInfo = true;
     else
         warning('.Date field is not formatted correctly')
         isSessionInfo = false;
     end
end

