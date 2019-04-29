function [positionTracking] = pos2behav(pos,varargin)
% [positionTracking] = pos2behav(pos,trackingType,varargin)
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
%    positionTracking       struct with the following fields
%                   .position      -reconstructed postions for each frame
%                            .x
%                            .y
%                            .z
%                            .units -unit of measurement ('meters', 'cm')
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
%                               .rotationType  -Default: 'quaternion', 
%                                               other option is euler
%                   .behaviorinfo
%                       .description
%                       .acquisitionsystem
%                       .substructnames
%                   .events        -Important time markers (behavioral scoring of trials, etc) 
%
%  Written by David Tingley, 2017

p=inputParser();
addParameter(p,'trials',[],@iscell);
addParameter(p,'mapping',[],@iscell);
addParameter(p,'map',[],@iscell);
addParameter(p,'behavType','',@isstr);
parse(p,varargin{:})

trials = p.Results.trials;
mapping = p.Results.mapping;
map = p.Results.map;
behavType = p.Results.behavType;

% determine if LED tracking or Motive tracking...
if size(pos,2) == 5
    trackingType = 'led';
    % check if timestamps are first or last
    [a b] = min(nanstd(diff(pos)));
    if b == 5
       pos = [pos(:,5),pos(:,1:4)];
       ts_col = 5;
       x_cols = [1 3];
       y_cols = [2 4];
    else
        ts_col = 1;
        x_cols = [2 4];
        y_cols = [3 5];
    end
elseif size(pos,2) == 11
    trackingType = 'optitrack';
else
    error('unrecognized pos size, should be Nx5 for LED or Nx11 for IR tracking')
end


switch trackingType
    case 'led'
        % check if timestamps are first or last column....
        [~,ts] = min(nanstd(diff(pos)));
        if ts == 1
            x_coords = nanmean(pos(:,[x_cols])');
            y_coords = nanmean(pos(:,[3 5])');
            [orientation rho] = cart2pol(pos(:,y_cols(1))-pos(:,y_cols(2)),pos(:,x_cols(2))-pos(:,x_cols(2)));
        else
            x_coords = nanmean(pos(:,[1 3])');
            y_coords = nanmean(pos(:,[x_cols])');
            [orientation rho] = cart2pol(pos(:,x_cols(1))-pos(:,x_cols(2)),pos(:,y_cols(1))-pos(:,y_cols(2)));
        end
        timestamps = pos(:,ts);
        
        positionTracking.position.x = x_coords';
        positionTracking.position.y = y_coords';
        positionTracking.position.z = [];
        positionTracking.timestamps = timestamps;
        positionTracking.samplingRate = 1000 ./ nanmean(diff(positionTracking.timestamps))./1000;
        positionTracking.position.units = 'pixels';
        positionTracking.orientation.yaw = orientation;
        positionTracking.orientation.pitch = [];
        positionTracking.orientation.roll = [];
        positionTracking.orientation.rotationType = 'euler';
        positionTracking.behaviorinfo.description = behavType;
        positionTracking.behaviorinfo.acquisitionsystem = trackingType;
        positionTracking.behaviorinfo.substructnames = {'position','orientation'};
        
        if ~isempty(trials)
            t=1;
            for i=1:length(trials)
                for j=1:length(trials{i})
                    positionTracking.events.trials{t}.x = nanmean(trials{i}{j}(:,x_cols)')';
                    positionTracking.events.trials{t}.y = nanmean(trials{i}{j}(:,y_cols)')';
                    positionTracking.events.trials{t}.z = [];
                    positionTracking.events.trials{t}.orientation.yaw = cart2pol(trials{i}{j}(:,x_cols(1))-...
                        trials{i}{j}(:,x_cols(2)),trials{i}{j}(:,y_cols(1))-trials{i}{j}(:,y_cols(2)));
                    positionTracking.events.trials{t}.orientation.pitch = [];
                    positionTracking.events.trials{t}.orientation.roll = [];
                    positionTracking.events.trials{t}.timestamps = trials{i}{j}(:,ts_col);
                    positionTracking.events.trials{t}.mapping = mapping{i}{j}(:,5);
                    map_c.x = nanmean(map{i}(:,x_cols)')';
                    map_c.y = nanmean(map{i}(:,y_cols)')';
                    map_c.z = [];
                    positionTracking.events.map{i} = map_c;
                    positionTracking.events.trialConditions(t) = i;
                    
                    positionTracking.events.trialIntervals(t,:) = [trials{i}{j}(1,5);trials{i}{j}(end,5)];
                    t = 1+t;
                end
            end
        end
                
    case 'optitrack'
        
        positionTracking.position.x = pos(:,8);
        positionTracking.position.y = pos(:,10);
        positionTracking.position.z = pos(:,9);
        positionTracking.timestamps = pos(:,1);
        positionTracking.samplingRate = 1000 ./ nanmean(diff(positionTracking.timestamps))./1000;
        if nanstd(pos(:,7)) > 10  % determine unit of measure
            positionTracking.position.units = 'meters';
        else
            positionTracking.position.units = 'mm';
        end
        positionTracking.orientation.x = pos(:,4);
        positionTracking.orientation.y = pos(:,5);
        positionTracking.orientation.z = pos(:,6);
        positionTracking.orientation.w = pos(:,7);
        positionTracking.orientation.rotationType = 'quaternion';
        positionTracking.behaviorinfo.errorPerMarker = pos(:,11);
        positionTracking.behaviorinfo.description = behavType;
        positionTracking.behaviorinfo.acquisitionsystem = trackingType;
        positionTracking.behaviorinfo.substructnames = {'position','orientation'};
        
         if ~isempty(trials)
            t=1;
            for i=1:length(trials)
                for j=1:length(trials{i})
                    positionTracking.events.trials{t}.x = trials{i}{j}(:,8);
                    positionTracking.events.trials{t}.y = trials{i}{j}(:,10);
                    positionTracking.events.trials{t}.z = trials{i}{j}(:,9);
                    positionTracking.events.trials{t}.orientation.x = trials{i}{j}(:,4);
                    positionTracking.events.trials{t}.orientation.y = trials{i}{j}(:,5);
                    positionTracking.events.trials{t}.orientation.z = trials{i}{j}(:,6);
                    positionTracking.events.trials{t}.orientation.w = trials{i}{j}(:,7);
                    positionTracking.events.trials{t}.errorPerMarker = trials{i}{j}(:,11);
                    positionTracking.orientation.rotationType = 'quaternion';
                    positionTracking.events.trials{t}.timestamps = trials{i}{j}(:,1);
                    positionTracking.events.trials{t}.mapping = mapping{i}{j}(:,13);
                    map_c.x = map{i}(:,8);
                    map_c.y = map{i}(:,10);
                    map_c.z = map{i}(:,9);
                    positionTracking.events.map{i} = map_c;
                    positionTracking.events.trialConditions(t) = i;
                    positionTracking.events.trialIntervals(t,:) = [trials{i}{j}(1,1);trials{i}{j}(end,1)];
                    t = 1+t;
                end
            end
        end
        
    otherwise
        error('unrecognized tracking type, have you added a new tracking method?')        
end





