function [behavior] = pos2behav(pos,trackingType,varargin)
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
        
        behavior.position.x = x_coords';
        behavior.position.y = y_coords';
        behavior.position.z = [];
        behavior.timestamps = timestamps;
        behavior.samplingRate = 1000 ./ mean(diff(behavior.timestamps))./1000;
        behavior.units = 'pixels';
        behavior.orientation.yaw = orientation;
        behavior.orientation.pitch = [];
        behavior.orientation.roll = [];
        behavior.rotationType = 'euler';
        behavior.description = behavType;
        behavior.trackingType = trackingType;
        
        if ~isempty(trials)
            t=1;
            for i=1:length(trials)
                for j=1:length(trials{i})
                    if size(trials{i}{j},2) == 5
                    behavior.events.trials{t}.x = nanmean(trials{i}{j}(:,x_cols)')';
                    behavior.events.trials{t}.y = nanmean(trials{i}{j}(:,y_cols)')';
                    behavior.events.trials{t}.z = [];
                    behavior.events.trials{t}.orientation.yaw = cart2pol(trials{i}{j}(:,x_cols(1))-...
                        trials{i}{j}(:,x_cols(2)),trials{i}{j}(:,y_cols(1))-trials{i}{j}(:,y_cols(2)));
                    behavior.events.trials{t}.orientation.pitch = [];
                    behavior.events.trials{t}.orientation.roll = [];
                    behavior.events.trials{t}.timestamps = trials{i}{j}(:,ts_col);
                    behavior.events.trials{t}.mapping = mapping{i}{j}(:,5);
                    map_c.x = nanmean(map{i}(:,x_cols)')';
                    map_c.y = nanmean(map{i}(:,y_cols)')';
                    map_c.z = [];
                    behavior.events.map{i} = map_c;
                    behavior.events.trialConditions(t) = i;
                    
                    behavior.events.trialIntervals(t,:) = [trials{i}{j}(1,5);trials{i}{j}(end,5)];
                    else
                        
                    behavior.events.trialIntervals(t,:) = [trials{i}{j}(1,1);trials{i}{j}(end,1)];
                    end
                    t = 1+t;
                end
            end
        end
                
    case 'optitrack'
        
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
        behavior.description = behavType;
        behavior.events = [];
        
    otherwise
        error('unrecognized tracking type, have you added a new tracking method?')        
end





