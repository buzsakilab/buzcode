datafolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/';
%Make List of all Recordings in data folder
recordings = dir(datafolder);
recordings(strncmpi('.',{recordings.name},1)) = [];   %Remove entries starting with .,~
recordings(strncmpi('~',{recordings.name},1)) = []; 

numrecs = length(recordings);
%%
auto = zeros(numrecs,1);
manS = zeros(numrecs,1);
manE = zeros(numrecs,1);
WS = zeros(numrecs,1);
WSE = zeros(numrecs,1);
autoint = zeros(numrecs,3);
manSint = zeros(numrecs,3);
manEint = zeros(numrecs,3);
WSint = zeros(numrecs,3);
WSEint = zeros(numrecs,3);
%%
for r = 1:numrecs;

display(['Recording ',num2str(r),' of ',num2str(numrecs)]);

recname = recordings(r).name;


load([datafolder,recname,'/',recname,'_Intervals.mat'])
oldints = intervals([1,3,5]);


load([datafolder,recname,'/',recname,'_StateIntervals.mat'])
WSstateints = {StateIntervals.Wake,StateIntervals.SWS,StateIntervals.REM};

load([datafolder,recname,'/',recname,'_StateIDM.mat'])
stateintsM = stateintervals;
episodeintsM = episodeintervals;

load([datafolder,recname,'/',recname,'_WSRestrictedIntervals.mat'])
WSepisodeints = {WakeInts,SWSEpisodeInts,REMEpisodeInts};

load([datafolder,recname,'/',recname,'_StateIDA.mat'])
stateintsA = stateintervals;

%%


[auto(r),autoint(r,:)] = CompareIntervalSets(stateintsA,oldints,'nulltimes',0);
[manS(r),manSint(r,:)] = CompareIntervalSets(stateintsM,oldints,'nulltimes',0);
[manE(r),manEint(r,:)] = CompareIntervalSets(episodeintsM,oldints,'nulltimes',0);
[WS(r),WSint(r,:)] = CompareIntervalSets(WSstateints,oldints,'nulltimes',0);
[WSE(r),WSEint(r,:)] = CompareIntervalSets(WSepisodeints,oldints,'nulltimes',0);

end



%%
statecolors = {'k','b','r'}
figure
    subplot(2,3,1)
        hold on
        boxplot([auto,autoint])
        plot(ones(size(auto)),auto,'.','color',0.7*[1 1 1],'MarkerSize',15)
        for ii = 1:3
            plot(1+ii*ones(size(auto)),autoint(:,ii),'.','color',statecolors{ii},'MarkerSize',20)
        end
        title('Automatic Scoring')
        ylabel('Similarity to Previous Scoring')
        xlim([0 5])
        set(gca,'XTick',[1:4])
        set(gca,'XTickLabel',{'All','W','N','R'})
        ylim([0 1])
        
    subplot(2,3,2)
        hold on
        boxplot([manS,manSint])
        plot(ones(size(manS)),manS,'.','color',0.7*[1 1 1],'MarkerSize',15)
        for ii = 1:3
            plot(1+ii*ones(size(auto)),manSint(:,ii),'.','color',statecolors{ii},'MarkerSize',15)
        end
        title('After Manual Check')
        ylabel('Similarity to Previous Scoring')
        xlim([0 5])
        set(gca,'XTick',[1:4])
        set(gca,'XTickLabel',{'All','W','N','R'})
        ylim([0 1])
    subplot(2,3,4.5)
        hold on
        boxplot([manE,manEint])
        plot(ones(size(manE)),manE,'.','color',0.7*[1 1 1],'MarkerSize',15)
        for ii = 1:3
            plot(1+ii*ones(size(auto)),manEint(:,ii),'.','color',statecolors{ii},'MarkerSize',15)
        end
        title('After Manual Check - Episodes')
        ylabel('Similarity to Previous Scoring')
        xlim([0 5])
        set(gca,'XTick',[1:4])
        set(gca,'XTickLabel',{'All','W','N','R'})
        ylim([0 1])
    subplot(2,3,3)
        hold on
        boxplot([WS,WSint])
        plot(ones(size(WS)),WS,'.','color',0.7*[1 1 1],'MarkerSize',15)
        for ii = 1:3
            plot(1+ii*ones(size(auto)),WSint(:,ii),'.','color',statecolors{ii},'MarkerSize',15)
        end
        title('Wake-Sleep for Analysis')
        ylabel('Similarity to Previous Scoring')
        xlim([0 5])
        set(gca,'XTick',[1:4])
        set(gca,'XTickLabel',{'All','W','N','R'})
        ylim([0 1])
    subplot(2,3,5.5)
        hold on
        boxplot([WSE,WSEint])
        plot(ones(size(WSE)),WSE,'.','color',0.7*[1 1 1],'MarkerSize',15)
        for ii = 1:3
            plot(1+ii*ones(size(auto)),WSEint(:,ii),'.','color',statecolors{ii},'MarkerSize',15)
        end
        title('Wake-Sleep for Analysis - Episodes')
        ylabel('Similarity to Previous Scoring')
        xlim([0 5])
        set(gca,'XTick',[1:4])
        set(gca,'XTickLabel',{'All','W','N','R'})
        ylim([0 1])
        
saveas(gcf,[figfolder,'oldvsnew'],'jpeg')