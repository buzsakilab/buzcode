%% POS conversion

posfile = dir('*pos');
pos = importdata(posfile.name);

% determine if timestamps are first or last column
[a b] = min(nanstd(diff(pos))); % find the column with smallest variability..

behav.timestamps = pos(:,b);
pos(:,b) = []; % remove timestamps from pos mat

if size(pos,2) > 5  % if optitrack
    columns = [7 9 8 3 5 4 6 10 1 2];
    pos = pos(:,columns);
    behav.position.x = pos(:,2);
    behav.position.y = pos(:,3);
    behav.position.z = pos(:,4);
    behav.orientation.rx = pos(:,5);
    behav.orientation.ry = pos(:,6);
    behav.orientation.rz = pos(:,7);
    behav.orientation.rw = pos(:,8);
    behav.timestamps = pos(:,1);
    behav.errorPerMarker = pos(:,9);
    behav.frameCount = pos(:,10);

elseif size(pos,2) < 5  % if LED tracking
    behav.position.x = nanmean(pos(:,[2 4]));
    behav.position.y = nanmean(pos(:,[1 3]));    
    dx = pos(:,3) - pos(:,5);
    dy = pos(:,2) - pos(:,4);
	ang = atan2(dy,dx)-angOffset;
	ang = mod(ang,2*pi);
    behav.orientation.z = ang; 
    warning('come up with a better head dir calculation...') 
end


%% LFP conversion



%% SPIKE conversion

% times

spktimes.timestamps = ;


save('spktimes.cellinfo.mat')

% waveforms


save('waveforms.cellinfo.mat')

% features

save('features.cellinfo.mat')

% metadata

save('metadata.cellinfo.mat')

%% METADATA conversion



%% EVENT conversion



%% move old files to FMAT_format folder



