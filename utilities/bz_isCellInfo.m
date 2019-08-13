function [iscellinfo] = bz_isCellInfo(cellinfo)
% USAGE
% [] = bz_isCellInfo(spikes)
% 
% INPUT
%       cellinfo  - struct with the following fields
%                   .UID
%                   .sessionName
%                   .region
%
% OUTPUT
%      logical true if struct meets cellinfo criteria, false if otherwise
%
% written by david tingley, 2017

if isfield(cellinfo,'UID') && isfield(cellinfo,'sessionName') && isfield(cellinfo,'region') % check that fields exist
     if ischar(cellinfo.sessionName) && isvector(cellinfo.UID) && isvector(cellinfo.region)
         iscellinfo = true;
     else
         warning('one of the required fields for a cellinfo type is not formatted correctly')
         iscellinfo = false;
     end
else
    warning('one of the required fields for a cellinfo type does not exist')
    iscellinfo = false;
end
