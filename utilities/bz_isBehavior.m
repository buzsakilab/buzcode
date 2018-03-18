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


if isfield(behavior,'behaviorinfo') && isfield(behavior,'timestamps') && isfield(behavior,'samplingRate')% check that fields exist
     if isstruct(behavior.behaviorinfo) && isvector(behavior.timestamps) && isnumeric(behavior.samplingRate)
         % check that sub fields exist
         if isfield(behavior.behaviorinfo,'description') && isfield(behavior.behaviorinfo,'acquisitionsystem') && isfield(behavior.behaviorinfo,'substructnames') 
            isBehavior = true;
         else
            warning('behavior struct seems to be missing something in .behaviorinfo')
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
