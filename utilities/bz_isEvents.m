<<<<<<< HEAD
function [isEvents] = bz_isEvents(events)
=======
function [isEvent] = bz_isEvents(event)
>>>>>>> dev
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


<<<<<<< HEAD
if isfield(events,'detectorinfo') && isfield(events,'timestamps') % check that fields exist
     if isstruct(events.detectorinfo) && isvector(events.timestamps) 
         % check that sub fields exist
         if isfield(events.detectorinfo,'detectorname') && isfield(events.detectorinfo,'detectionparms') ...
                 && isfield(events.detectorinfo,'detectioninterval') && isfield(events.detectorinfo,'detectiondate')
            isEvents = true;
         else
             isEvents = false;
=======
if isfield(event,'detectorinfo') && isfield(event,'timestamps') % check that fields exist
     if isstruct(event.detectorinfo) && isvector(event.timestamps) 
         % check that sub fields exist
         if isfield(event.detectorinfo,'detectorname') && isfield(event.detectorinfo,'detectionparms') ...
                 && isfield(event.detectorinfo,'detectioninterval') && isfield(event.detectorinfo,'detectiondate')
            isEvent = true;
         else
             isEvent = false;
>>>>>>> dev
         end
         
     else
         warning('one of the required fields for an event type is not formatted correctly')
<<<<<<< HEAD
         isEvents = false;
     end
else
    warning('one of the required fields for an event type does not exist')
    isEvents = false;
=======
         isEvent = false;
     end
else
    warning('one of the required fields for an event type does not exist')
    isEvent = false;
>>>>>>> dev
end
