function [isBehavior] = bz_isBehavior(behavior)
% USAGE
% [isBehavior] = bz_isBehavior(behavior)
% 
% INPUT
%       behavior   - struct with the following fields
%                   .timestamps
%                   .position
%                       .detectorname
%                       .detectionparms
%                       .detectioninterval
%                       .detectiondate
%
% OUTPUT
%      logical true if struct meets behavior criteria, false if otherwise
%
% written by david tingley, 2017


if isfield(behavior,'position') && isfield(behavior,'units') && isfield(behavior,'timestamps') && isfield(behavior,'samplingRate')% check that fields exist
     if isstruct(behavior.position) && isvector(behavior.timestamps) && ischar(behavior.units) && isnumeric(behavior.samplingRate)
         % check that sub fields exist
         if isfield(behavior.position,'x') && isfield(behavior.position,'y') 
            isBehavior = true;
         else
            isBehavior = false;
         end
         
     else
         warning('one of the required fields for an behavior type is not formatted correctly')
         isBehavior = false;
     end
else
    warning('one of the required fields for an behavior type does not exist')
    isBehavior = false;
end
