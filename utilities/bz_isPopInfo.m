function [isPopInfo] = bz_isPopInfo(PopInfo)
% USAGE
% [] = bz_isPopInfo(spikes)
% 
% INPUT
%       PopInfo  - struct with the following fields
%                   .UID
%                   .sessionName
%                   .region
%
% OUTPUT
%      logical true if struct meets PopInfo criteria, false if otherwise
%
% written by david tingley, 2017


if isfield(PopInfo,'UID') && isfield(PopInfo,'sessionName') && isfield(PopInfo,'region') % check that fields exist
     if ischar(PopInfo.sessionName) && isvector(PopInfo.UID) && isvector(PopInfo.region)
         isPopInfo = true;
     else
         warning('one of the required fields for a PopInfo type is not formatted correctly')
         isPopInfo = false;
     end
else
    warning('one of the required fields for a PopInfo type does not exist')
    isPopInfo = false;
end
