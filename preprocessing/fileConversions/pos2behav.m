function [] = pos2behav(pos,varargin)
%pos2behav - convert old .pos files/variables to buzcode 'behavior' struct.
%
%  USAGE
%    [behav] = pos2beahv(varargin)
%
% INPUTS
%    pos            Nx5 position data (columns 2 and 4 are X, 3 and 5 are
%                   Y, and 1 is timestamps; there are two X/Y variables for two LED's)
%
%    OR
%
%    pos            Nx11 position data (optitrack)
%                     1 - Recording time (seconds, relative to recording start)
%                     2 - Frame count
%                     3 - Frame time (seconds, relative to tracking start)
%                     4 - X rotation
%                     5 - Y rotation
%                     6 - Z rotation
%                     7 - W rotation (quaternion rotation)
%                     8 - X position
%                     9 - Y position
%                     10- Z position
%                     11- Error per marker
%                     NOTE: Optitrack systems use Y-coordinate system,
%                     Y = height, X and Z = horizontal plane
%
%   trackingType    string indentifier of the tracking type, currently
%                   accepted inputs are 'led', 'optitrack', and ...
%
%  OUTPUT
%    behavior       struct with the following fields
%                   .position      -reconstructed postions for each frame
%                            .x
%                            .y
%                            .z
%                   .units         -unit of measurement ('meters', 'cm')
%                   .timestamps    -Nx1 vector of timestamps (relative to
%                                   recording) that match the x,y,z positions
%                   .samplingRate
%                   .orientation
%                               .x
%                               .y
%                               .z
%                               .w
%                               OR
%                               .pitch
%                               .yaw
%                               .roll
%                   .rotationType  -Default: 'quaternion', other option is
%                                   euler
%                   .description
%                   .events        -Important time markers (behavioral scoring of trials, etc) 
%
%  Written by David Tingley, 2017


% determine if LED tracking or Motive tracking...
if size(pos,2) == 5
    trackingType = 'led';
elseif size(pos,2) == 11
    trackingType = 'optitrack';
else
    error('unrecognized pos size, should be Nx5 for LED or Nx11 for IR tracking')
end


switch trackingType
    case trackingType == 'led'
        % check if timestamps are first or last column....
        [~,ts] = min(nanstd(diff(pos)));
        if ts == 1
            x_coords = nanmean(pos(:,[2 4])');
            y_coords = nanmean(pos(:,[3 5])');
            [orientation rho] = cart2pol(pos(:,3)-pos(:,5),pos(:,2)-pos(:,4));
        else
            x_coords = nanmean(pos(:,[1 3])');
            y_coords = nanmean(pos(:,[2 4])');
            [orientation rho] = cart2pol(pos(:,2)-pos(:,4),pos(:,1)-pos(:,3));
        end
        timestamps = pos(:,ts);
        
        behavior.position.x = x_coords;
        behavior.position.y = y_coords;
        behavior.position.z = [];
        behavior.timestamps = timestamps;
        behavior.samplingRate = 1000 ./ mean(diff(behavior.timestamps))./1000;
        behavior.units = 'pixels';
        behavior.orientation.yaw = orientation;
        behavior.orientation.pitch = [];
        behavior.orientation.roll = [];
        behavior.rotationType = 'euler';
        behavior.description = '';
        behavior.events = [];
        
    case trackingType == 'optitrack'
        
        behavior.position.x = pos(:,8);
        behavior.position.y = pos(:,10);
        behavior.position.z = pos(:,9);
        behavior.timestamps = pos(:,1);
        behavior.samplingRate = 1000 ./ mean(diff(behavior.timestamps))./1000;
        if nanstd(pos(:,7)) > 10  % determine unit of measure
            behavior.units = 'meters';
        else
            behavior.units = 'mm';
        end
        behavior.orientation.x = pos(:,4);
        behavior.orientation.y = pos(:,5);
        behavior.orientation.z = pos(:,6);
        behavior.orientation.w = pos(:,7);
        behavior.rotationType = 'quaternion';
        behavior.errorPerMarker = pos(:,11);
        behavior.description = '';
        behavior.events = [];
        
    otherwise
        error('unrecognized tracking type, have you added a new tracking method?')        
end





