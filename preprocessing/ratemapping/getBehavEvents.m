function [x y pos trials map mapping] = getBehavEvents(pos,bins)
% this function takes unlabeled position tracking data (LED or Optitrack) and,
% through a semi-automated process, scores individual trials where an animal 
% traverses through the environment.  Specifically, this code works well with the
% 2AFC T-maze, lap alternation on the circular maze, or traversals on a linear track

% assumptions about the 'pos' variable need to be worked out... (introduced by LED vs Optitrack)


% David Tingley, 2017 

if pos(1,1) == 0 & size(pos,2) < 6
    %% old code for LED tracking
    pos = [pos(:,2:5),pos(:,1)];    
    % assumes colomns 1/2 and 3/4 match
    scatter(pos(:,1),pos(:,2),'.r')
    hold on
    scatter(pos(:,3),pos(:,4),'.b')
%     axis([-1 1 -1 1])
    pos(intersect(find(pos(:,1)==0),find(pos(:,3)==0)),:) = nan;
    pos(intersect(find(pos(:,2)==0),find(pos(:,4)==0)),:) = nan;


    p(:,1) = nanmean(pos(:,[1 3])');
    f = find(isnan(p(:,1)));
    p(:,2) = nanmean(pos(:,[2 4])');
    
    velThresh = 8; 
    distThresh = 50; 
    x_coord = 1;
    y_coord = 2;
    
    vel = nansum(abs(diff(p)'))';
    vel(vel>100) = 0;
    vel = fastrms(vel,4);
    subplot(2,1,1)
    hist(vel(~isnan(p(2:end,1))),100);
    axis([ 0 mean(vel)+8*std(vel) 0 max(hist(vel(~isnan(p(2:end,1))),100))])
    subplot(2,1,2)
    plot(vel(~isnan(p(2:end,1)))',1:sum(~isnan(p(2:end,1)))')
    axis([ 0 mean(vel)+8*std(vel) 0 sum(~isnan(p(2:end,1)))])
%     thresh = min(-vel);
%     while max(diff((lo))) > 500
%         [pl lo]=findpeaks(-vel,'MINPEAKHEIGHT',thresh);
%         thresh = thresh - min(-vel)./100 ;
%     end
    [x y]=ginput();
    velThresh=x; %7e-5; % change this to first min peak in vel histogram
    close all
    scatter(p(:,1),p(:,2),1,'.k')
    distThresh =  mean(max(p)-min(p))./2;
end
if size(pos,2)>6
    %% new code for IR Motive tracking
%     if std(pos(:,y_coord)) > 50
%        pos(:,x_coord) = pos(:,x_coord)./1000;
%        pos(:,9) = pos(:,9)./1000;
%        pos(:,y_coord) = pos(:,y_coord)./1000; 
%     end
    x_coord = 8;
    y_coord = 10;
    p = pos(:,[x_coord y_coord]);
    if nanstd(pos(:,x_coord)) > 100
        p(:,1) = p(:,1)./1000;
        
        p(:,2) = p(:,2)./1000;
    end
%     p(isnan(p(:,1)),:)=[];
    m = min(p(:));
    p = p-m;
%     p(:,1) = fastrms(p(:,1),12);
%     p(:,2) = fastrms(p(:,2),12);
    p(:,1) = smoothts(p(:,1),'b',40);
    p(:,2) = smoothts(p(:,2),'b',40);

    vel = nansum(abs(diff(p)'))';
    vel(vel>100) = 0;
    vel = fastrms(fastrms(vel,4),30);
    subplot(2,1,1)
    hist(vel(~isnan(p(2:end,1))),1000);
    axis([ 0 mean(vel)+8*std(vel) 0 max(hist(vel(~isnan(p(2:end,1))),1000))])
    subplot(2,1,2)
    plot(vel(~isnan(p(2:end,1)))',1:sum(~isnan(p(2:end,1)))')
    axis([ 0 mean(vel)+8*std(vel) 0 sum(~isnan(p(2:end,1)))])
%     thresh = min(-vel);
%     while max(diff((lo))) > 500
%         [pl lo]=findpeaks(-vel,'MINPEAKHEIGHT',thresh);
%         thresh = thresh - min(-vel)./100 ;
%     end
    [x y]=ginput();
    velThresh=x; %7e-5; % change this to first min peak in vel histogram
    close all
    scatter(p(:,1),p(:,2),1,'.k')
    distThresh =  mean(max(p)-min(p))./20;

%     axis([-1 1 -1 1])
end
dbstop if error
 


disp('pick start/stop locations')
[x,y] = ginput();
x = x(2:end);
y = y(2:end); % laptop hack....
ff = find(isnan(p(:,2)));
% p=pos(:,1:2);

for i=1:length(x)
    dists(i,:) = abs(p(:,1)-x(i))+abs(p(:,2)-y(i));
    [pks{i} locs{i}]=findpeaks(-dists(i,:),'MinPeakHeight',-distThresh);
    locs{i}=locs{i}';
    for j=1:length(x)
    trials{i,j,1}=[];
    trials{i,j,2}=[];
    end
end
[a a locsmat] = spikes2sorted(locs);
f = find(locsmat(:,1)<120);
locsmat(f,:)=[];

centerx = mean(x);
centery = mean(y);

lastStop = 0;
lastStart = 0;
clf();
cc=1;

for l=1:length(locsmat)
    next=1;
    while l+next+12 < length(locsmat) & locsmat(l,1) ~= locsmat(l+next,1)
        % check that the velocity doesn't drop below threshold
        if any(vel(locsmat(l,1):locsmat(l+next,1))<velThresh)
            break
        end
        % check that the run is long enough to be at least one lap 
%         if (sum(sum(abs(diff(p(locsmat(l,1):locsmat(l+next,:),:)))))) > distThresh
            
            % check that the velocity drops below threshold within a
            % certain number of frames
%             if any(vel(locsmat(l,1)-120:locsmat(l,1)+120)<velThresh) & any(vel(locsmat(l+next,1)-120:locsmat(l+next,1)+120)<velThresh)
            startPos=pos(locsmat(l,1),:);
            stopPos=pos(locsmat(l+next,1),:);
            if lastStop < startPos(1) % cut off four seconds for some 'wiggle' room
%                 [a b]=min(abs(nanmean(startPos([1 3]))-x)+abs(nanmean(startPos([2 4]))-y));
%                 [aa c]=min(abs(nanmean(stopPos([1 3]))-x)+abs(nanmean(stopPos([2 4]))-y)); 
                [a b]=min(abs((startPos(x_coord))-x)+abs((startPos(y_coord))-y));
                [aa c]=min(abs((stopPos(x_coord))-x)+abs((stopPos(y_coord))-y));  % set if for old vs new tracking system here
                    subplot(2,2,1)
                    hold off
                    plot(p(:,1),p(:,2),'.k')
                    hold on
                    scatter(p(locsmat(l,1)-80:locsmat(l+next,1)+80,1),p(locsmat(l,1)-80:locsmat(l+next,1)+80,2),300,'.b')
                    scatter(p(locsmat(l,1),1),p(locsmat(l,1),2),300,'.g')
                    scatter(p(locsmat(l+next,1),1),p(locsmat(l+next,1),2),300,'.r')
                    subplot(2,2,2)
                    plot(vel(locsmat(l,1)-80:locsmat(l+next,1)+80))
                    hold on
                    line([20 20],[0 velThresh*2],'color','g')
                    line([locsmat(l+next,1)-locsmat(l,1) locsmat(l+next,1)-locsmat(l,1)],[0 velThresh*2],'color','r')
                    hold off
                    
                a = input('all good [w/a/s/d/y/n]: ','s');
%                 a=temp{cc}; 
%                cc=1+cc;
%                 length(temp)
                while ~strcmp(a,'') & ~strcmp(a,'n') & ~strcmp(a,'y')
                    if strcmp(a,'d')
                        next=1+next;
                    end
                    if strcmp(a,'a')
                        l=l-1;
                        next=next+1;
                    end
                    if strcmp(a,'w')
                        l=l+1;
                        next=next-1;
                    end
                    if strcmp(a,'s')
                        next=next-1;
                    end
                    subplot(2,2,1)
                    hold off
                    plot(p(:,1),p(:,2),'.k')
                    hold on
                    scatter(p(locsmat(l,1)-80:locsmat(l+next,1)+80,1),p(locsmat(l,1)-80:locsmat(l+next,1)+80,2),300,'.b')
                    scatter(p(locsmat(l,1),1),p(locsmat(l,1),2),300,'.g')
                    scatter(p(locsmat(l+next,1),1),p(locsmat(l+next,1),2),300,'.r')
                    subplot(2,2,2)
                    plot(vel(locsmat(l,1)-80:locsmat(l+next,1)+80))
                    hold on
                    line([20 20],[0 velThresh*2],'color','g')
                    line([locsmat(l+next,1)-locsmat(l,1) locsmat(l+next,1)-locsmat(l,1)],[0 velThresh*2],'color','r')
                    hold off
                    a = input('all good [w/a/s/d/y/n]: ','s');
%                     a=temp{cc}; 
%                     cc=cc+1;
                end
                if strcmp(a,'y') | strcmp(a,'') % add good trial
                    a = p(locsmat(l+next,1)-100:locsmat(l+next,1)-20,1);
                    bb = p(locsmat(l+next,1)-100:locsmat(l+next,1)-20,2);
                    a(isnan(a))=[];
                    bb(isnan(bb))=[];
                pp=polyfit([a],...
                        [bb],1);    
              trials{b,c,heaviside(pp(1))+1}{end+1}=pos(locsmat(l,1):locsmat(l+next,1),:);
                lastStart = startPos(1); 
                lastStop = stopPos(1);
                end
                if strcmp(a,'n')
                    lastStart = startPos(1); 
                    lastStop = stopPos(1);
                end
            end
%             end
%         end
        next=next+1;
    end   
end


trials = reshape(trials,length(x)*length(x)*2,1); 
c=1;
for i=1:length(trials)
    for j=1:length(trials{i})
        trials_unsorted{c} = trials{i}{j};
        trials_new{c} = trials{i}{j};
        c=1+c;
    end
end
%% merge trials to the same length
for i=1:length(trials_unsorted)
       d = pdist(trials_unsorted{i}(:,[x_coord y_coord]),'euclidean');
       dd = squareform(d);
       dd(dd==0)=nan;
   while length(trials_unsorted{i})>bins

       [minVal]=min(dd(:));
       [row col] = find(minVal == dd);
       row = row(end);
       col = col(end);
       % check that row and col are adjacent 
       if abs((row-col))<=1
           trials_unsorted{i} = [trials_unsorted{i}(1:row-1,:);...
                                mean(trials_unsorted{i}(row:col,:));...
                                trials_unsorted{i}(col+1:end,:)];
           d = pdist(trials_unsorted{i}(:,[x_coord y_coord]),'euclidean');
           dd = squareform(d);
           dd(dd==0)=nan;
       elseif abs((row-col))>1
           dd(row,col)=nan;
       end
   end
end
%% TODO, add automatic removal of trials that are too short
% cluster similarity
for i=1:length(trials_unsorted)
    for j =1:length(trials_unsorted)
        a = trials_unsorted{i}(:,[x_coord y_coord]);
        b = trials_unsorted{j}(:,[x_coord y_coord]);
        temp = corrcoef(reshape(a,bins*2,1),reshape(b,bins*2,1));
        cc(i,j)=temp(2);
    end
end

%% get number of trial types / behavioral condition manually...
numConditions = str2num(input('How many conditions were there: ','s'));
exitC=0;
while exitC == 0
    clust = kmeans(cc,numConditions);
%     clust = GMM(cc,numConditions);
    clear trials
    for i=1:numConditions
        f = find(clust==i);
        for j=1:length(f)
            trials{i}{j}=trials_new{f(j)};
        end
    end
    for i=1:numConditions
        subplot(4,5,i)
        hold off
        for t = 1:length(trials{i})
            scatter(trials{i}{t}(:,x_coord),trials{i}{t}(:,y_coord),'.k')
            hold on
            scatter(trials{i}{t}(1,x_coord),trials{i}{t}(1,y_coord),'.g')
            scatter(trials{i}{t}(end,x_coord),trials{i}{t}(end,y_coord),'.r')
        end
    end
    pause(.01);
    adjust = input('merge[m]/split[s] ','s');
    if strcmp(adjust,'s')
        numConditions = numConditions+1;
    end
    while strcmp(adjust,'m')
        one = str2num(input('condition #1:','s'));
        two = str2num(input('condition #2:','s'));
        for temp = 1:length(trials{two})
           trials{one}{length(trials{one})+1}=trials{two}{temp};
        end
        trials{two}=[];
        trials = trials(~cellfun('isempty',trials));  
        numConditions = numConditions-1;
        clf
        for i=1:numConditions
            subplot(4,5,i)
            hold off
            for t = 1:length(trials{i})
                scatter(trials{i}{t}(:,x_coord),trials{i}{t}(:,y_coord),'.k')
                hold on
                scatter(trials{i}{t}(1,x_coord),trials{i}{t}(1,y_coord),'.g')
                scatter(trials{i}{t}(end,x_coord),trials{i}{t}(end,y_coord),'.r')
            end
        end
       adjust = input('merge[m]/split[s] ','s');
    end
    if ~strcmp(adjust,'s') && ~strcmp(adjust,'m')
       exitC = 1; 
       clf();
    end
end
% check that clustering was correct 

%% normalize positions to template
c=1;
for tt = 1:length(trials)
%             if d == 5  %% only for rec 21st
%                 bins = 100;
%             elseif d == 4
%                 bins = 100;
%             else
%                 bins = 80;
%             end
    
    hold off
    map{tt}=[];
    t_conc=[];
    for t = 1:length(trials{tt})
        t_conc = [trials{tt}{t}(:,:),20*(trials{tt}{t}(:,1)-trials{tt}{t}(1,1))];
    if length(t_conc)>=bins
        while length(t_conc)>bins+1
        di = pdist(t_conc);
        s = squareform(di);
        s(find(eye(size(s))))=nan;
        [a b] = min(s(:));
        [coords blah] = find(s==a);
        t_conc(coords(1),:) = (t_conc(coords(1),:)+t_conc(coords(2),:))./2;
        t_conc(coords(2),:) = [];
        % debug
%         scatter(t_conc(:,1),t_conc(:,2))
%         pause(.01)
        end
    t_conc_all(t,:,:) = t_conc;
    end
    end
    if length(trials{tt})>0
    map{tt} = squeeze(median(t_conc_all(:,:,:),1));

    if ~isempty(trials{tt})
    subplot(8,4,c)
    c=1+c
    for t = 1:length(trials{tt})
%     scatter(map{tt}(:,1),map{tt}(:,2),'.')
    scatter(trials{tt}{t}(:,1),trials{tt}{t}(:,2),'.k')
    hold on
    axis([0 550 0 550])
    end
    for t = 1:length(trials{tt})
%     scatter(map{tt}(:,1),map{tt}(:,2),'.')
    scatter(trials{tt}{t}(1,1),trials{tt}{t}(1,2),'.g')
    scatter(trials{tt}{t}(end,1),trials{tt}{t}(end,2),'.r')
    hold on
    
    axis([0 550 0 550])
    end
    end
    end
    clear t_conc_all
%     for t = 1:length(trials{tt})
%         scatter((trials{tt}{t}(:,1)),(trials{tt}{t}(:,2)),'.r')
%         hold on
%         en
%         scatter((trials{tt}{t}(:,3)),(trials{tt}{t}(:,4)),'.b')
%         scatter((trials{tt}{t}(1,3)),(trials{tt}{t}(1,4)),'.g')
%     end
%     if ~isempty(trials{tt})
%     g = ginput();
%     ggx=[];ggy=[];
%     for t =1:length(g)-1whps
%         dd = abs(g(t,1)-g(t+1,1)) + abs(g(t,2)-g(t+1,2));
%         ggx = [ggx,linspace(g(t,1),g(t+1,1),round(dd))];
%     
%         ggy = [ggy,linspace(g(t,2),g(t+1,2),round(dd))];
%     end
%     g = [ggx;ggy]';
%     
%     step = (max(g(:,2))-min(g(:,2)))./bins;
%     step2 = (max(g(:,2))-min(g(:,2)))./length(g(:,2));
%     ww = smooth(interp1(min(g(:,2)):step2+.0001:max(g(:,2)),g(:,2),min(g(:,2)):step:max(g(:,2))+.01,'linear','extrap'),3);
%     step = (max(g(:,1))-min(g(:,1)))./bins;
%     step2 = (max(g(:,1))-min(g(:,1)))./length(g(:,1));
%     w = smooth(interp1(min(g(:,1)):step2+.0001:max(g(:,1)),g(:,1),min(g(:,1)):step:max(g(:,1))+.01,'linear','extrap'),3);
%     plot(w,ww,'k')
    disp('finding mapping...')
%     map{tt} = [w,ww];
    for t =1:length(trials{tt})  % all trial types (rotations)
        for p = 1:length(trials{tt}{t})
            [a b] = min(nansum(abs([trials{tt}{t}(p,1)-map{tt}(:,1),...
                trials{tt}{t}(p,x_coord)-map{tt}(:,x_coord),...
                trials{tt}{t}(p,9)-map{tt}(:,9),...
                trials{tt}{t}(p,y_coord)-map{tt}(:,y_coord),...
                (trials{tt}{t}(p,1)-trials{tt}{t}(1,1))*50-map{tt}(:,1),...  % penalty for time differences
                40*(p./length(trials{tt}{t})*length(map{tt}) - (1:length(map{tt})))'])'));     % penalty for order differences
            mapping{tt}{t}(p,:) = [map{tt}(b,1:end) b trials{tt}{t}(p,1)];
%             plot(nansum(abs([trials{tt}{t}(p,1)-map{tt}(:,1),...
%                 trials{tt}{t}(p,2)-map{tt}(:,2),...
%                 trials{tt}{t}(p,3)-map{tt}(:,3),...
%                 trials{tt}{t}(p,4)-map{tt}(:,4),...
%                 (trials{tt}{t}(p,5)-trials{tt}{t}(1,5))*20-map{tt}(:,5),...  % penalty for time differences
%                 50*(p./length(trials{tt}{t})*length(map{tt}) - (1:length(map{tt})))'])'))
%             pause
        end
    end
end

%% incorporate this below to reformat things and double check trial assignments...
% load behav
% 
c=1
clear m mm tt
for i=1:length(trials)
    if ~isempty(trials{i})
    tt{c}=trials{i};
    m{c}=map{i};
    mm{c}=mapping{i};
    
    tt{c}=tt{c}(~cellfun('isempty',tt{c}));
    mm{c}=mm{c}(~cellfun('isempty',mm{c}));
    c=1+c;
    end
end 
dbmap=m;
mapping=mm;
trials= tt;

for tt = 1:length(trials)
subplot(5,4,tt)
for t = 1:length(trials{tt})
%     scatter(map{tt}(:,1),map{tt}(:,2),'.')
scatter(trials{tt}{t}(:,x_coord),trials{tt}{t}(:,y_coord),'.k')
hold on

%     axis([0 550 0 550])
end
for t = 1:length(trials{tt})
%     scatter(map{tt}(:,1),map{tt}(:,2),'.')
scatter(trials{tt}{t}(1,x_coord),trials{tt}{t}(1,y_coord),'.g')
scatter(trials{tt}{t}(end,x_coord),trials{tt}{t}(end,y_coord),'.r')
hold on

%     axis([0 550 0 550])
end
title(tt);
end

%% update to return buzcode format
% behavior = pos2behav(pos,'optitrack','trials',trials,'mapping',mapping,'map',map,'behavType','central alternation');