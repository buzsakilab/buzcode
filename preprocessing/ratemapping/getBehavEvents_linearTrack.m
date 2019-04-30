function [x y pos trials map mapping] = getBehavEvents_linearTrack(pos,bins)
dbstop if error
scatter(pos(:,8),pos(:,10),'.')
% axis([-1 1 -1 1])
[x y]=ginput();
cx =  x(end-1:end);
y = y(end-1:end);
dist(:,1)=fastrms(abs(pos(:,8)-x(1))+abs(pos(:,10)-y(1)),120);
dist(:,2)=fastrms(abs(pos(:,8)-x(2))+abs(pos(:,10)-y(2)),120);

plot(dist,'.')
% axis([0 length(dist) 0 3])
[xx yy] = ginput();

[pks_one locs_one]=findpeaks(-dist(:,1),'MINPEAKHEIGHT',-yy);
[pks_two locs_two]=findpeaks(-dist(:,2),'MINPEAKHEIGHT',-yy);
c=1;
lastTrial = 0;
for i=1:length(locs_one)
    f = locs_one(i)-locs_two;
    ff = find(f>0);
    [a b]=min(locs_one(i)-locs_two(ff));
    if a < 500 & a > 60 & locs_two(ff(b)) > lastTrial + 60
       trials{c} = pos(locs_two(ff(b)):locs_one(i),:); 
       c=1+c;
       lastTrial = locs_two(ff(b));
    end
    f = locs_two-locs_one(i);
    ff = find(f>0);
    [a b]=min(locs_two(ff)-locs_one(i));
    if a < 500 & a > 60 & locs_one(i) > lastTrial + 60
       trials{c} = pos(locs_one(i):locs_two(ff(b)),:); 
       c=1+c;
       lastTrial = locs_one(i);
    end
end

% for i=1:length(locs_one)
%     two = locs_two(locs_two-locs_one(i)>0);
%     if ~isempty(two)
%     [a b]=min(abs(two-locs_one(i)));  
%     if two(b)-locs_one(i) < 24000 && two(b)-locs_one(i) > 120 && locs_one(i) > lastTrial-60
%         trials{c} = pos(locs_one(i):two(b),:);
%         c=1+c;
%         lastTrial = two(b);
%     end
%     end
%     two = locs_two(locs_one(i)-locs_two>0);
%     if ~isempty(two)
%     [a b]=min(abs(two-locs_one(i)));  
%     if locs_one(i)-two(b) < 24000 && locs_one(i)-two(b) > 120 && two(b) > lastTrial-60
%         trials{c} = pos(two(b):locs_one(i),:);
%         c=1+c;
%         lastTrial = locs_one(i);
%     end
%     end
% end

%% fine tune trials here
% for i=1:length(trials)
% clear t tt ttt tttt; t = trials{i};
% tt=t(:,8:10);
% d = diff(trials{i}(:,8:10));
% tt(find(d==0)+1)=nan;
% for j=1:3
% tttt(:,j) = smooth(fillmissing(tt(:,j),'pchip'),15);
% end
% subplot(3,1,1)
% plot((tttt))
% axis([0 length(tttt) -1 1])
% % hold on
% % plot(tttt(:,2));
% % hold off
% subplot(3,1,2)
% plot(diff(tttt))
% hold on
% plot(sum(abs(diff(tttt(:,[1 3])))'))
% hold off
% axis([0 length(tttt) -.015 .015])
% subplot(3,1,3)
% plot(abs(tttt(:,1)-tttt(:,3)))
% axis([0 length(tttt) 0 1.4])
% % scatter3(tttt(:,1),tttt(:,3),tttt(:,2),'.k')
% [xxx yyy]=ginput();
% if xxx(2) > length(trials{i})
%     xxx(2)=length(trials{i});
% end
% trials_new{i}=trials{i}(round(xxx(1)):round(xxx(2)),:)
% end


%% merge trials to the same length
trials_unsorted = trials;
for i=1:length(trials_unsorted)
       d = pdist(trials_unsorted{i}(:,[8 9 10]),'euclidean');
       dd = squareform(d);
       dd(dd==0)=nan;
       [mins ind] = min(dd(:));
           [row col] = ind2sub(size(dd),ind);
%            f = find(abs(row-col)==1);  % this worked but it was sloowww,
%            row = row(f);               % removed and replaced w ind2sub
%            col = col(f);
    while length(trials_unsorted{i})>bins
         [mins ind] = min(dd(:));
         [row col] = ind2sub(size(dd),ind);
       % check that row and col are adjacent 
       if ((row(1)-col(1)))== -1
           trials_unsorted{i} = [trials_unsorted{i}(1:row(1)-1,:);...
                                nanmean(trials_unsorted{i}(row(1):col(1),:));...
                                trials_unsorted{i}(col(1)+1:end,:)];          
           dd(row,:)=[];
           dd(:,col)=[];
       elseif ((row(1)-col(1)))== 1
           trials_unsorted{i} = [trials_unsorted{i}(1:col(1)-1,:);...
                                nanmean(trials_unsorted{i}(col(1):row(1),:));...
                                trials_unsorted{i}(row(1)+1:end,:)];
           dd(row,:)=[];
           dd(:,col)=[];
       elseif abs((row(1)-col(1)))>1 | abs((row(1)-col(1)))<1
           dd(row,col)=nan;
       end
    end
    i
end

% cluster similarity
for i=1:length(trials_unsorted)
    for j =1:length(trials_unsorted)
        for ax = 8:10
        temp(ax-7,:)=crosscorr(trials_unsorted{i}(:,[ax]),trials_unsorted{j}(:,[ax]),bins-1);
        end
        if sum(~isreal(temp(:)))>0
            disp([i j]);
        end
        a = trials_unsorted{i}(:,[8 9 10]);
        b = trials_unsorted{j}(:,[8 9 10]);
% a = [trials_unsorted{i}(:,[9]) trials_unsorted{i}(:,[10])-trials_unsorted{i}(:,[9])];
% b = [trials_unsorted{j}(:,[9]) trials_unsorted{j}(:,[10])-trials_unsorted{j}(:,[9])];
%         temp = corrcoef(reshape(a,bins*3,1),reshape(b,bins*2,1));
%           temp = crosscorr(reshape(a,bins*3,1),reshape(b,bins*,1),120);
%         cc(i,j)=temp(2);
          [cc(i,j) locs]=max(nanmean(temp));
    end
end
% remove NaN's?
n = isnan(cc(:,1));

%% 
numConditions = str2num(input('How many conditions were there: ','s'));
while ~isnumeric(numConditions) || isempty(numConditions)
    numConditions = str2num(input('How many conditions were there: ','s'));
end
exitC=0;
while exitC == 0
    clust = kmeans(cc,numConditions);
    clear trials
    for i=1:numConditions
        f = find(clust==i);
        for j=1:length(f)
            trials{i}{j}=trials_unsorted{f(j)};
        end
    end
    for i=1:numConditions
        subplot(4,5,i)
        hold off
        for t = 1:length(trials{i})
            scatter3(trials{i}{t}(:,8),trials{i}{t}(:,10),trials{i}{t}(:,9),'.k')
            hold on
            scatter3(trials{i}{t}(1,8),trials{i}{t}(1,10),trials{i}{t}(1,9),'.g')
            scatter3(trials{i}{t}(end,8),trials{i}{t}(end,10),trials{i}{t}(end,9),'.r')
        end
    end
    pause(.01);
    adjust = input('merge[m]/split[s] ','s');
    if strcmp(adjust,'s')
        numConditions = numConditions+1;
    end
    if strcmp(adjust,'m')
        on = str2num(input('condition #1:','s'));
        tw = str2num(input('condition #2','s'));
        for t = 1:length(trials{tw})
           trials{on}{length(trials{on})+1}=trials{tw}{t};
        end
        trials{tw}=[];
    end
    if ~strcmp(adjust,'s') && ~strcmp(adjust,'m')
       exitC = 1; 
       clf();
    end
end
% check that clustering was correct 

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
    scatter(trials{tt}{t}(:,1),trials{tt}{t}(:,2),'.k')
    hold on
    scatter(trials{tt}{t}(1,1),trials{tt}{t}(1,2),'.g')
    scatter(trials{tt}{t}(end,1),trials{tt}{t}(end,2),'.r')
    end
    end
    end
    clear t_conc_all
    for t =1:length(trials{tt})  % all trial types (rotations)
        for p = 1:length(trials{tt}{t})
            [a b] = min(nansum(abs([trials{tt}{t}(p,1)-map{tt}(:,1),...
                trials{tt}{t}(p,8)-map{tt}(:,8),...
                trials{tt}{t}(p,9)-map{tt}(:,9),...
                trials{tt}{t}(p,10)-map{tt}(:,10),...
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
scatter(trials{tt}{t}(:,8),trials{tt}{t}(:,10),'.k')
hold on

%     axis([0 550 0 550])
end
for t = 1:length(trials{tt})
%     scatter(map{tt}(:,1),map{tt}(:,2),'.')
scatter(trials{tt}{t}(1,8),trials{tt}{t}(1,10),'.g')
scatter(trials{tt}{t}(end,8),trials{tt}{t}(end,10),'.r')
hold on

%     axis([0 550 0 550])
end
title(tt);
end


end