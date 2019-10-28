function [ ISIstats ] = bz_ISIStats( spikes,varargin )
%ISIstats = bz_ISIStats(spikes,varargin) calculates the statistics 
%inter-spike intervals for the spiketimes in spikes.
%
%   INPUTS
%       spikes      Structure with spikes.times and spikes.UID
%
%       (options)
%       'ints'        A structure with intervals in which to calculate ISIs.
%                       states.stateNAME = [start stop]
%                       Will calculate ISIsstats separately for each state
%                       (Can also 'load' from SleepState.states.mat)
%       'cellclass'     Cell array of strings - label for each cell. 
%                       (Can also 'load' from CellClass.cellinfo.mat)
%       'savecellinfo'  logical (default: false) save a cellinfo file?
%       'basePath'
%       'figfolder'     a folder to save the figure in
%       'showfig'       logical (default: false) show the figure?
%       'forceRedetect' logical (default: false) to re-compute even if saved
%       'shuffleCV2'    logical (devault: false)
%       'numISIbins'    number of bins for ISI distribution (default: 60)
%       'numCV2bins'    number pf bins for CV2 distribution (default: 60)
%
%   OUTPUTS
%       ISIstats        cellinfo structure with ISI statistics
%           .summstats  summary statistics
%           .ISIhist    histograms of ISIs etc
%           .sorts      sorting indices
%           .allspikes  ISI/CV2 value for each spike 
%                       (ISI is PRECEDING interval for the spike at allspikes.times)
%
%DLevenstein 2018
%% Parse the inputs
defaultstates.ALL = [-Inf Inf];

% parse args
p = inputParser;
addParameter(p,'ints',defaultstates)
addParameter(p,'savecellinfo',false,@islogical)
addParameter(p,'basePath',pwd,@isstr)
addParameter(p,'figfolder',false)
addParameter(p,'showfig',false,@islogical);
addParameter(p,'cellclass',[]);
addParameter(p,'forceRedetect',false,@islogical);
addParameter(p,'shuffleCV2',false,@islogical);
addParameter(p,'numISIbins',60,@islogical);
addParameter(p,'numCV2bins',50,@islogical);


parse(p,varargin{:})
ints = p.Results.ints;
cellclass = p.Results.cellclass;
basePath = p.Results.basePath;
SAVECELLINFO = p.Results.savecellinfo;
figfolder = p.Results.figfolder;
SHOWFIG = p.Results.showfig;
forceRedetect = p.Results.forceRedetect;
SHUFFLECV2 = p.Results.shuffleCV2;
numISIbins = p.Results.numISIbins;
numCV2bins = p.Results.numCV2bins;


%% Load the stuff
baseName = bz_BasenameFromBasepath(basePath);
cellinfofilename = fullfile(basePath,[baseName,'.ISIStats.cellinfo.mat']);

if exist(cellinfofilename,'file') && ~forceRedetect
    ISIstats = bz_LoadCellinfo(basePath,'ISIStats');
    return
end
    
if strcmp(cellclass,'load')
    cellclass = bz_LoadCellinfo(basePath,'CellClass');
    cellclass = cellclass.label;
end

if strcmp(ints,'load')
    ints = bz_LoadStates(basePath,'SleepState');
    ints = ints.ints;
end


%% Get the States
ints.ALL = [-Inf Inf];  %Add ALL to the options
statenames = fieldnames(ints);
numstates = length(statenames);


%% ISI and CV2 statistics
numcells = length(spikes.times);

%Calculate ISI and CV2 for allspikes
allspikes.ISIs = cellfun(@diff,spikes.times,'UniformOutput',false);
allspikes.meanISI = cellfun(@(X) (X(1:end-1)+X(2:end))./2,allspikes.ISIs,'UniformOutput',false);
allspikes.CV2 = cellfun(@(X) 2.*abs(X(2:end)-X(1:end-1))./(X(2:end)+X(1:end-1)),allspikes.ISIs ,'UniformOutput',false);
%Make sure times line up
allspikes.times = cellfun(@(X) X(2:end-1),spikes.times,'UniformOutput',false);
allspikes.ISIs = cellfun(@(X) X(1:end-1),allspikes.ISIs,'UniformOutput',false);
%%
for ss = 1:numstates
%ss=1;

%Find which spikes are during state of interest
statespikes = cellfun(@(X) InIntervals(X,ints.(statenames{ss})),...
    allspikes.times,'UniformOutput',false);
CV2 = cellfun(@(X,Y) X(Y(2:end-1)),allspikes.CV2,statespikes,'Uniformoutput',false);
ISIs = cellfun(@(X,Y) X(Y(2:end-1)),allspikes.ISIs,statespikes,'Uniformoutput',false);
normISIs = cellfun(@(X) X./mean(X),ISIs,'Uniformoutput',false);


%Summary Statistics
summstats.(statenames{ss}).meanISI = cellfun(@(X) mean(X),ISIs);
summstats.(statenames{ss}).meanrate = 1./summstats.(statenames{ss}).meanISI;
summstats.(statenames{ss}).ISICV = cellfun(@(X) std(X)./mean(X),ISIs);
summstats.(statenames{ss}).meanCV2 = cellfun(@(X) mean(X),CV2);



%Account for no spikes
summstats.(statenames{ss}).meanrate(isnan(summstats.(statenames{ss}).meanrate))=0;


if SHUFFLECV2
    % CV2 for shuffled ISIs
    numshuffle = 100;
    for sh = 1:numshuffle
        ISIs_shuffle = cellfun(@(X) shuffle(X),ISIs,'UniformOutput',false);
        CV2_shuffle = cellfun(@(X) 2.*abs(X(2:end)-X(1:end-1))./(X(2:end)+X(1:end-1)),...
            ISIs_shuffle ,'UniformOutput',false);
        meanshuffle(sh,:) = cellfun(@(X) mean(X),CV2_shuffle);
    end
    summstats.(statenames{ss}).shufflemeanCV2 = mean(meanshuffle);
    summstats.(statenames{ss}).shufflestdCV2 = std(meanshuffle);

end
%%
%Set up all the bins and matrices
ISIhist.linbins = linspace(0,10,numISIbins);
ISIhist.logbins = linspace(log10(0.001),log10(200),numISIbins);
ISIhist.(statenames{ss}).lin = zeros(numcells,numISIbins);
ISIhist.(statenames{ss}).log = zeros(numcells,numISIbins);
ISIhist.(statenames{ss}).return = zeros(numISIbins,numISIbins,numcells);
ISIhist.(statenames{ss}).meannorm = zeros(numcells,numISIbins);

CV2hist.bins = linspace(0,2,numCV2bins+1);
CV2hist.bins = CV2hist.bins(1:end-1)+0.5.*diff(CV2hist.bins([1 2]));
CV2hist.(statenames{ss}) = zeros(numcells,numCV2bins);

Jointhist.(statenames{ss}).log = zeros(numcells,numISIbins,numCV2bins);
Jointhist.(statenames{ss}).norm = zeros(numcells,numISIbins,numCV2bins);



%Calculate all the histograms: ISI, log(ISI), 1/ISI, log(1/ISI)
for cc = 1:numcells
    numspks(cc) = length(ISIs{cc});
    spikethresh = 100; %need 100 spikes to show in mean plot
    enoughspikes(cc) = numspks(cc)>spikethresh;
    
    if numspks(cc)<1
        continue
    end
    %Calculate ISI histograms
    ISIhist.(statenames{ss}).lin(cc,:) = hist(ISIs{cc},ISIhist.linbins);
    ISIhist.(statenames{ss}).log(cc,:) = hist(log10(ISIs{cc}),ISIhist.logbins);
    ISIhist.(statenames{ss}).loginv(cc,:) = hist(log10(1./ISIs{cc}),ISIhist.logbins);
    ISIhist.(statenames{ss}).meannorm(cc,:) = hist(log10(normISIs{cc}),ISIhist.logbins);
    
    CV2hist.(statenames{ss})(cc,:) = hist(CV2{cc},CV2hist.bins);
    Jointhist.(statenames{ss}).log(cc,:,:) = hist3(...
        [log10([ISIs{cc};ISIs{cc}(2:end)]),[CV2{cc};CV2{cc}(1:end-1)]],...
        {ISIhist.logbins,CV2hist.bins});
    Jointhist.(statenames{ss}).norm(cc,:,:) = hist3(...
        [log10([normISIs{cc};normISIs{cc}(2:end)]),[CV2{cc};CV2{cc}(1:end-1)]],...
        {ISIhist.logbins,CV2hist.bins});
    
    
    %Normalize histograms to number of spikes
    ISIhist.(statenames{ss}).lin(cc,:) = ISIhist.(statenames{ss}).lin(cc,:)./numspks(cc);
    ISIhist.(statenames{ss}).log(cc,:) = ISIhist.(statenames{ss}).log(cc,:)./numspks(cc);
    ISIhist.(statenames{ss}).loginv(cc,:) = ISIhist.(statenames{ss}).loginv(cc,:)./numspks(cc);
    ISIhist.(statenames{ss}).meannorm(cc,:) = ISIhist.(statenames{ss}).meannorm(cc,:)./numspks(cc);
    
    CV2hist.(statenames{ss})(cc,:) = CV2hist.(statenames{ss})(cc,:)./numspks(cc);
    
    Jointhist.(statenames{ss}).log(cc,:,:) = Jointhist.(statenames{ss}).log(cc,:,:)./numspks(cc);
    Jointhist.(statenames{ss}).norm(cc,:,:) = Jointhist.(statenames{ss}).norm(cc,:,:)./numspks(cc);

    %Calculate Return maps
    if numspks(cc)>1
    ISIhist.(statenames{ss}).return(:,:,cc) = hist3(log10([ISIs{cc}(1:end-1) ISIs{cc}(2:end)]),{ISIhist.logbins,ISIhist.logbins});
    end
    ISIhist.(statenames{ss}).return(:,:,cc) = ISIhist.(statenames{ss}).return(:,:,cc)./numspks(cc);
  
end

summstats.(statenames{ss}).numspikes = numspks;

%Sortings
[~,sorts.(statenames{ss}).rate]=sort(summstats.(statenames{ss}).meanrate);
[~,sorts.(statenames{ss}).ISICV]=sort(summstats.(statenames{ss}).ISICV);
[~,sorts.(statenames{ss}).CV2]=sort(summstats.(statenames{ss}).meanCV2);

%Make the cell-type specific sortings and average distributions
if ~isempty(cellclass)
    %Check for empty cell class entries
    noclass = cellfun(@isempty,cellclass);
    sorts.numclassycells = sum(~noclass);
    %cellclass(noclass)={'none'};
    classnames = unique(cellclass(~noclass));
    numclasses = length(classnames);
    for cl = 1:numclasses
        inclasscells{cl} = strcmp(classnames{cl},cellclass);
        
        %Mean distributions
        meandists.(statenames{ss}).(classnames{cl}).ISIdist = squeeze(nanmean(ISIhist.(statenames{ss}).log(inclasscells{cl}&enoughspikes,:),1));
        meandists.(statenames{ss}).(classnames{cl}).CV2dist = squeeze(nanmean(CV2hist.(statenames{ss})(inclasscells{cl}&enoughspikes,:),1));
        meandists.(statenames{ss}).(classnames{cl}).Jointdist = squeeze(nanmean(Jointhist.(statenames{ss}).log(inclasscells{cl}&enoughspikes,:,:),1));
        
        %Sorts
        sorttypes = {'rate','ISICV','CV2'};
        for tt = 1:length(sorttypes)
        sorts.(statenames{ss}).([sorttypes{tt},classnames{cl}]) = ...
            intersect(sorts.(statenames{ss}).(sorttypes{tt}),find(inclasscells{cl}),'stable');
        
        if cl==1
            sorts.(statenames{ss}).([sorttypes{tt},'byclass'])=[];
        end
        sorts.(statenames{ss}).([sorttypes{tt},'byclass']) = ...
            [sorts.(statenames{ss}).([sorttypes{tt},'byclass']) sorts.(statenames{ss}).([sorttypes{tt},classnames{cl}])];
        end
            
    end  
else
    numclasses = 1; cl=1;
    classnames{1} = 'ALL';
    sorts.numclassycells = numcells;
    inclasscells{1} = true(1,numcells);
    sorts.(statenames{ss}).ratebyclass = sorts.(statenames{ss}).rate;
    sorts.(statenames{ss}).ISICVbyclass = sorts.(statenames{ss}).ISICV;
    sorts.(statenames{ss}).CV2byclass = sorts.(statenames{ss}).CV2;
    
        %Mean distributions
        meandists.(statenames{ss}).(classnames{cl}).ISIdist = squeeze(nanmean(ISIhist.(statenames{ss}).log(inclasscells{cl}&enoughspikes,:),1));
        meandists.(statenames{ss}).(classnames{cl}).CV2dist = squeeze(nanmean(CV2hist.(statenames{ss})(inclasscells{cl},:)&enoughspikes,1));
        meandists.(statenames{ss}).(classnames{cl}).Jointdist = squeeze(nanmean(Jointhist.(statenames{ss}).log(inclasscells{cl}&enoughspikes,:,:),1));
end



%% ACG - get its own analysis
%[ccg,t] = CCG(statespiketimes,[],<options>)

%%
if SHOWFIG | figfolder
figure
    subplot(2,2,1)
        for cl = 1:numclasses
            plot(log10(summstats.(statenames{ss}).meanrate(inclasscells{cl})),...
                log2(summstats.(statenames{ss}).ISICV(inclasscells{cl})),'.','markersize',11)
            hold on
        end
        plot(get(gca,'xlim'),log2([1 1]),'k')
        LogScale('x',10);LogScale('y',2);
        xlabel('Mean Rate (Hz)');ylabel('ISI CV')
        title(statenames{ss})
        box off
        
    subplot(2,2,2)
    for cl = 1:numclasses
        
        if SHUFFLECV2
            plot([log10(summstats.(statenames{ss}).meanrate(inclasscells{cl}));...
                log10(summstats.(statenames{ss}).meanrate(inclasscells{cl}))],...
                [summstats.(statenames{ss}).shufflemeanCV2(inclasscells{cl});...
                summstats.(statenames{ss}).meanCV2(inclasscells{cl})],...
                'color',0.8.*[1 1 1],'linewidth',0.25)
            hold on
            plot([log10(summstats.(statenames{ss}).meanrate(inclasscells{cl}));...
                log10(summstats.(statenames{ss}).meanrate(inclasscells{cl}))],...
                [summstats.(statenames{ss}).shufflemeanCV2(inclasscells{cl})-summstats.(statenames{ss}).shufflestdCV2(inclasscells{cl});...
                summstats.(statenames{ss}).shufflemeanCV2(inclasscells{cl})+summstats.(statenames{ss}).shufflestdCV2(inclasscells{cl})],...
                'color',0.8.*[1 1 1],'linewidth',2.5)
        end
        
        plot(log10(summstats.(statenames{ss}).meanrate(inclasscells{cl})),...
            (summstats.(statenames{ss}).meanCV2(inclasscells{cl})),'.','markersize',11)
        hold on
        
        

        
    end
        plot(get(gca,'xlim'),[1 1],'k')
        LogScale('x',10);
        xlabel('Mean Rate (Hz)');ylabel('ISI <CV2>')
        title(statenames{ss})
        box off

    subplot(2,3,4)
        imagesc((ISIhist.logbins),[1 sorts.numclassycells],...
            ISIhist.(statenames{ss}).log(sorts.(statenames{ss}).ratebyclass,:))
        hold on
        plot(log10(1./(summstats.(statenames{ss}).meanrate(sorts.(statenames{ss}).ratebyclass))),[1:sorts.numclassycells],'k.','LineWidth',2)
        plot(ISIhist.logbins([1 end]),sum(inclasscells{1}).*[1 1]+0.5,'r')
        LogScale('x',10)
        xlabel('ISI (s)')
        xlim(ISIhist.logbins([1 end]))
        %colorbar
      %  legend('1/Mean Firing Rate (s)','location','southeast')
        ylabel('Cell (Sorted by FR, Type)')
        %legend('1/Mean Firing Rate (s)','location','southeast')
        caxis([0 0.1])
        title('ISI Distribution (Log Scale)')
        
        
    subplot(2,3,5)
        imagesc((CV2hist.bins),[1 sorts.numclassycells],...
            CV2hist.(statenames{ss})(sorts.(statenames{ss}).ratebyclass,:))
        hold on
        plot((summstats.(statenames{ss}).meanCV2(sorts.(statenames{ss}).ratebyclass)),[1:sorts.numclassycells],'k.','LineWidth',2)
        plot(CV2hist.bins([1 end]),sum(inclasscells{1}).*[1 1]+0.5,'r')
        %LogScale('x',10)
        xlabel('CV2')
        xlim(CV2hist.bins([1 end]))
        %colorbar
        xlim([0 2])
      %  legend('1/Mean Firing Rate (s)','location','southeast')
        ylabel('Cell (Sorted by FR, Type)')
        %legend('1/Mean Firing Rate (s)','location','southeast')
        %caxis([0 0.1])
        title('CV2 Distribution')
        
	for cl = 1:numclasses
        subplot((numclasses.*2),3,cl.*3+(numclasses.*3))
            imagesc(ISIhist.logbins,CV2hist.bins,meandists.(statenames{ss}).(classnames{cl}).Jointdist')
            hold on
            plot(ISIhist.logbins,meandists.(statenames{ss}).(classnames{cl}).ISIdist.*20,'k')
            plot(meandists.(statenames{ss}).(classnames{cl}).CV2dist.*20-3,CV2hist.bins,'k')
            axis xy
            LogScale('x',10)
            ylabel({(classnames{cl}),'CV2'});
            xlabel('ISI (s)');
            
            bottomrow = meandists.(statenames{ss}).(classnames{cl}).Jointdist(:,1);
            caxis([0 max([bottomrow;0])])
	end
        
%     subplot(2,3,5)
%         imagesc((ISIhist.logbins),[1 sorts.numclassycells],...
%             ISIhist.(statenames{ss}).log(sorts.(statenames{ss}).ISICVbyclass,:))
%         hold on
%         plot(log10(1./(summstats.(statenames{ss}).meanrate(sorts.(statenames{ss}).ISICVbyclass))),[1:sorts.numclassycells],'k.','LineWidth',2)
%         plot(ISIhist.logbins([1 end]),sum(inclasscells{1}).*[1 1]+0.5,'r')
%         LogScale('x',10)
%         xlabel('ISI (s)')
%         xlim(ISIhist.logbins([1 end]))
%         %colorbar
%       %  legend('1/Mean Firing Rate (s)','location','southeast')
%         ylabel('Cell (Sorted by CV, Type)')
%         %legend('1/Mean Firing Rate (s)','location','southeast')
%         caxis([0 0.1])
%         title('ISI Distribution (Log Scale)')
%         
%     subplot(2,3,6)
%         imagesc((ISIhist.logbins),[1 sorts.numclassycells],...
%             ISIhist.(statenames{ss}).log(sorts.(statenames{ss}).CV2byclass,:))
%         hold on
%         plot(log10(1./(summstats.(statenames{ss}).meanrate(sorts.(statenames{ss}).CV2byclass))),[1:sorts.numclassycells],'k.','LineWidth',2)
%         plot(ISIhist.logbins([1 end]),sum(inclasscells{1}).*[1 1]+0.5,'r')
%         LogScale('x',10)
%         xlabel('ISI (s)')
%         xlim(ISIhist.logbins([1 end]))
%         %colorbar
%       %  legend('1/Mean Firing Rate (s)','location','southeast')
%         ylabel('Cell (Sorted by CV2, Type)')
%         %legend('1/Mean Firing Rate (s)','location','southeast')
%         caxis([0 0.1])
%         title('ISI Distribution (Log Scale)')
        
%     subplot(2,2,3)
%         plot(log2(summstats.(statenames{ss}).meanrate(CellClass.pE)),...
%             log2(summstats.(statenames{ss}).meanCV2(CellClass.pE)),'k.')
%         hold on
%         plot(log2(summstats.(statenames{ss}).meanrate(CellClass.pI)),...
%             log2(summstats.(statenames{ss}).meanCV2(CellClass.pI)),'r.')
%         LogScale('xy',2)
%         xlabel('Mean Rate (Hz)');ylabel('Mean CV2')
%         title(statenames{ss})
%         box off


if figfolder
    NiceSave(['ISIstats_',(statenames{ss})],figfolder,baseName);
end


%exneurons %top/bottom 25 %ile rate/ISICV
%%
% exwindur = 4; %s
% STATEtimepoints = Restrict(lfp.timestamps,double(SleepState.ints.(statenames{ss})));
% samplewin = STATEtimepoints(randi(length(STATEtimepoints))) + [0 exwindur];
% %%
% figure
% bz_MultiLFPPlot( lfp,'spikes',spikes,'timewin',samplewin,...
%     'sortmetric',summstats.(statenames{ss}).meanCV2,...
%     'cellgroups',{CellClass.pI,CellClass.pE})


end

ISIstats.summstats = summstats;
ISIstats.ISIhist = ISIhist;
ISIstats.CV2hist = CV2hist;
ISIstats.Jointhist = Jointhist;
ISIstats.sorts = sorts;
ISIstats.UID = spikes.UID;
ISIstats.allspikes = allspikes;


if SAVECELLINFO
    save(cellinfofilename,'ISIstats')
end
    
end

