function [isEvents] = bz_isEvents(events)
% USAGE
% [isEvent] = bz_isEvents(event)
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


if isfield(events,'detectorinfo') && isfield(events,'timestamps') % check that fields exist
     if isstruct(events.detectorinfo) && isnumeric(events.timestamps) 
         % check that sub fields exist
         if isfield(events.detectorinfo,'detectorname') && isfield(events.detectorinfo,'detectionparms') ...
                 && isfield(events.detectorinfo,'detectionintervals') && isfield(events.detectorinfo,'detectiondate')
            isEvents = true;
         else
             warning('one of the required fields for an event type is not formatted correctly')
             isEvents = false;
         end
         
     else
         warning('one of the required fields for an event type is not formatted correctly')
         isEvents = false;
     end
else
    warning('one of the required fields for an event type does not exist')
    isEvents = false;
end
