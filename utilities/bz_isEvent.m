function [isEvent] = bz_isEvent(event)
% USAGE
% [isEvent] = bz_isEvent(event)
% 
% INPUT
%       event   - struct with the following fields
%                   .timestamps
%                   .detectorinfo
%                       .detectorname
%                       .detectionparms
%                       .detectioninterval
%                       .detectiondate
%
% OUTPUT
%      logical true if struct meets event criteria, false if otherwise
%
% written by david tingley, 2017


if isfield(event,'detectorinfo') && isfield(event,'timestamps') % check that fields exist
     if isstruct(event.detectorinfo) && isvector(event.timestamps) 
         % check that sub fields exist
         if isfield(event.detectorinfo,'detectorname') && isfield(event.detectorinfo,'detectionparms') ...
                 && isfield(event.detectorinfo,'detectioninterval') && isfield(event.detectorinfo,'detectiondate')
            isEvent = true;
         else
             isEvent = false;
         end
         
     else
         warning('one of the required fields for an event type is not formatted correctly')
         isEvent = false;
     end
else
    warning('one of the required fields for an event type does not exist')
    isEvent = false;
end
