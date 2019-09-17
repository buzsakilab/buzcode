function [behavior] = bz_getJumpBehav(pos)


if nanstd(pos(:,8)) < 10
    pos(:,8:10) = pos(:,8:10)*1000;    
end
p = pos(:,[8 10]);
behavior = pos2behav(pos,'behavType','jump');

p(:,1) = smoothts(p(:,1),'b',40);
p(:,2) = smoothts(p(:,2),'b',40);

% subplot(2,2,1)
% scatter(pos(:,10),pos(:,8),'.')
% [x y] = ginput();

subplot(2,1,1)
scatter(pos(:,8)-pos(:,10),pos(:,9),'.')
axis([min(pos(:,8)-pos(:,10)) max(pos(:,8)-pos(:,10)) min(pos(:,9)) max(pos(:,9))])
[x y] = ginput();

for i=1:length(x)
    dists_x(i,:) = abs(pos(:,8)-pos(:,10)-x(i));
    dists_y(i,:) = abs(pos(:,9)-y(i));
    
    [pks_x{i} locs_x{i}] = findpeaks(-dists_x(i,:),'MinPeakHeight',-500);
    [pks_y{i} locs_y{i}] = findpeaks(-dists_y(i,:),'MinPeakHeight',-50);
    for j=1:length(x)
    trials{i,j,1}=[];
    trials{i,j,2}=[];
    end
    jump_ts{i} = 0;
    for j=length(pks_x{i}):-1:1
        for k=length(pks_y{i}):-1:1
            if abs(pks_x{i}(j)) < 150 & abs(pks_y{i}(k)) < 20
               if abs(locs_x{i}(j) - locs_y{i}(k)) < 10 
                   if abs(mean([locs_x{i}(j) locs_y{i}(k)]) - jump_ts{i}(end)) > 150
                    jump_ts{i} = [jump_ts{i};round(mean([locs_x{i}(j) locs_y{i}(k)]))];
                   end
               end
            end            
        end
    end
    jump_ts{i} = jump_ts{i}(2:end);
end

vel = nansum(abs(diff(p)'))';
vel = fastrms(fastrms(vel,40),80);
count = 1;
for k=1:length(x)
    for i = 1:length(jump_ts{k})
        if jump_ts{k}(i)+250 < length(vel)
        subplot(2,2,1)
        scatter(pos(:,8),pos(:,10),'.')
        hold on
        scatter(pos(jump_ts{k}(i)-40:jump_ts{k}(i)+40,8),pos(jump_ts{k}(i)-40:jump_ts{k}(i)+40,10),'.r')
        hold off
        subplot(2,2,2)
        scatter(pos(:,8),pos(:,9),'.')
        hold on
        scatter(pos(jump_ts{k}(i)-40:jump_ts{k}(i)+40,8),pos(jump_ts{k}(i)-40:jump_ts{k}(i)+40,9),'.r')
        axis([min(pos(:,8)-pos(:,10)) max(pos(:,8)-pos(:,10)) min(pos(:,9)) max(pos(:,9))])
        hold off
        subplot(2,2,3)
        plot(vel(jump_ts{k}(i)-250:jump_ts{k}(i)+250))
        subplot(2,2,4)
        plot(diff(vel(jump_ts{k}(i)-250:jump_ts{k}(i)+250)))
        hold on
        plot(medfilt1(diff(diff(vel(jump_ts{k}(i)-250:jump_ts{k}(i)+250)))*10,5))
        hold off
        
        [maxVel maxVelLoc] = max((vel(jump_ts{k}(i)-250:jump_ts{k}(i)+250)));
        [maxAcc maxAccLoc] = max(diff(vel(jump_ts{k}(i)-250:jump_ts{k}(i)+250)));
        [maxImp maxImpLoc] = max(medfilt1(diff(diff(vel(jump_ts{k}(i)-250:jump_ts{k}(i)+250))),5));
        if maxVelLoc ~= 501
        maxVelLoc = maxVelLoc + jump_ts{k}(i)-250;
        maxAccLoc = maxAccLoc + jump_ts{k}(i)-250;
        maxImpLoc = maxImpLoc + jump_ts{k}(i)-250;
        if maxVelLoc > maxAccLoc & maxAccLoc > maxImpLoc
            pause; s{k}(i) = 'y';
%         s{k}(i) = input('keep: [Y/n]','s');
        if strcmp(s{k}(i),'y')
            behavior.events.trials{count}.x =pos(maxImpLoc-400:maxImpLoc+400,8);
            behavior.events.trials{count}.y =pos(maxImpLoc-400:maxImpLoc+400,10);
            behavior.events.trials{count}.z =pos(maxImpLoc-400:maxImpLoc+400,9);
            behavior.events.trials{count}.mapping = 1:801;
            behavior.events.trials{count}.timestamps = pos(maxImpLoc-400:maxImpLoc+400,1);
            behavior.events.trials{count}.errorPerMarker = pos(maxImpLoc-400:maxImpLoc+400,11);

            
%             [r1 r2 r3] = quat2angle([pos(maxImpLoc-400:maxImpLoc+400,4);pos(maxImpLoc-400:maxImpLoc+400,5)...
%                 pos(maxImpLoc-400:maxImpLoc+400,6);pos(maxImpLoc-400:maxImpLoc+400,7)])
            behavior.events.trials{count}.orientation.x = pos(maxImpLoc-400:maxImpLoc+400,4);
            behavior.events.trials{count}.orientation.y = pos(maxImpLoc-400:maxImpLoc+400,5);
            behavior.events.trials{count}.orientation.z = pos(maxImpLoc-400:maxImpLoc+400,6);
            behavior.events.trials{count}.orientation.w = pos(maxImpLoc-400:maxImpLoc+400,7);

            if sum(diff(behavior.events.trials{count}.x)) > 0 & sum(diff(behavior.events.trials{count}.x)) > 0
               mod = 1;
            else
               mod = 2;
            end
            behavior.events.trialConditions(count) = k * mod;
            behavior.events.trialIntervals(count,:) = pos([maxImpLoc-400 maxImpLoc+400],1);
            count = 1+count;
        end
        end
        end
        end
    end
    behavior.events.conditionType{k} = 'jump';
end













