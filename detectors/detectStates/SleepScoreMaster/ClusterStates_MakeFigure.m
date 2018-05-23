function ClusterStates_MakeFigure(SleepState,basePath) 

%Check needed here if Sleep state has the right stuff
v2struct(SleepState.detectorinfo.detectionparms.SleepScoreMetrics)
v2struct(SleepState.detectorinfo.StatePlotMaterials)
v2struct(SleepState.detectorinfo.detectionparms.SleepScoreMetrics.histsandthreshs)

%
states = {SleepState.ints.WAKEstate,SleepState.ints.NREMstate,SleepState.ints.REMstate};
%%
%Figure locations
figloc = [fullfile(basePath,'StateScoreFigures'),'/'];
if ~exist(figloc,'dir')
    mkdir(figloc)
end

%% Figure
[zFFTspec,mu,sig] = zscore(log10(swFFTspec)');
[~,mu_th,sig_th] = zscore(log10(thFFTspec)');

 viewwin  =[t_clus(1) t_clus(end)];
 %viewwin  =[32000 34000];
%viewwin=[9000 11000];
clusterfig = figure('visible','off');
	subplot(8,1,[1:2])
        imagesc(t_clus,log2(swFFTfreqs),log10(swFFTspec))
        axis xy
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        caxis([3.5 6.5])
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        xlim(viewwin)
        colorbar('east')
        ylim([log2(swFFTfreqs(1)) log2(swFFTfreqs(end))+0.2])
        set(gca,'XTickLabel',{})
        ylabel({'swLFP','f (Hz)'})
        title([recordingname,': State Scoring Results']);
	subplot(8,1,3)
        imagesc(t_clus,log2(thFFTfreqs),log10(thFFTspec))
        axis xy
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        %caxis([3.5 6.5])
        caxis([min(mu_th)-2.5*max(sig_th) max(mu_th)+2.5*max(sig_th)])
        xlim(viewwin)
        %colorbar('east')
        ylim([log2(thFFTfreqs(1)) log2(thFFTfreqs(end))+0.2])
        ylabel({'thLFP','f (Hz)'})
        set(gca,'XTickLabel',{})
        
    subplot(8,1,4)
        %plot(t_clus,-IDX,'LineWidth',2)
        hold on
        plot(states{1}',-1*ones(size(states{1}))','k','LineWidth',8)
        plot(states{2}',-2*ones(size(states{2}))','b','LineWidth',8)
        plot(states{3}',-3*ones(size(states{3}))','r','LineWidth',8)
        xlim(viewwin)
        ylim([-4 0])
        set(gca,'YTick',[-3:-1])
        set(gca,'YTickLabel',{'REM','SWS','Wake/MA'})
        set(gca,'XTickLabel',{})
        
   	subplot(6,1,4)
        hold on
        plot(t_clus,broadbandSlowWave,'k')
        %plot(synchtimes',thresh*ones(size(synchtimes))','r')
        ylabel('SW')
        xlim([t_clus(1) t_clus(end)])
        xlim(viewwin)
        set(gca,'XTickLabel',{})
        
   	subplot(6,1,5)
        hold on
        plot(t_clus,thratio,'k')
        %plot(synchtimes',thresh*ones(size(synchtimes))','r')
        ylabel('Theta')
        xlim([t_clus(1) t_clus(end)])
        xlim(viewwin)
        set(gca,'XTickLabel',{})
        
   	subplot(6,1,6)
        hold on
        plot(t_EMG,EMG,'k')
        %plot(synchtimes',thresh*ones(size(synchtimes))','r')
        ylabel('EMG')
        xlim([t_clus(1) t_clus(end)])
        xlim(viewwin)
        xlabel('t (s)')
        
	saveas(clusterfig,[figloc,recordingname,'_SSResults'],'jpeg')
        

        
%% Figure: Split REM/Arousal  
IDX = bz_INTtoIDX(states,'length',max(t_clus));
IDX(1:t_clus(1)-1)=[];
NREMtimes = (broadbandSlowWave >swthresh);

figure
	subplot(3,2,1)
        hold on
        bar(swhistbins(swhistbins>swthresh),swhist(swhistbins>swthresh),'FaceColor','b','barwidth',0.9,'linewidth',1)
        bar(swhistbins(swhistbins<=swthresh),swhist(swhistbins<=swthresh),'FaceColor',0.9*[1 1 1],'barwidth',0.9,'linewidth',1)
        plot([swthresh swthresh],[0 max(swhist)],'r','LineWidth',1)
        xlabel('Broadband Slow Wave')
        title('Step 1: Broadband for NREM')
        

	subplot(3,2,3)
        hold on
        bar(EMGhistbins(EMGhistbins>EMGthresh),EMGhist(EMGhistbins>EMGthresh),'FaceColor','k','barwidth',0.9,'linewidth',1)
        bar(EMGhistbins(EMGhistbins<=EMGthresh),EMGhist(EMGhistbins<=EMGthresh),'FaceColor',0.9*[1 1 1],'barwidth',0.9,'linewidth',1)
        plot([EMGthresh EMGthresh],[0 max(EMGhist)],'r','LineWidth',1)
        xlabel('EMG')
        title('Step 2: EMG for Muscle Tone')
	subplot(3,2,5)
        hold on
        bar(THhistbins(THhistbins>=THthresh),THhist(THhistbins>=THthresh),'FaceColor','r','barwidth',0.9,'linewidth',1)
        bar(THhistbins(THhistbins<THthresh),THhist(THhistbins<THthresh),'FaceColor','k','barwidth',0.9,'linewidth',1)
        plot([THthresh THthresh],[0 max(THhist)],'r','LineWidth',1)
        xlabel('Theta')
        title('Step 3: Theta for REM')
        
        
    subplot(2,2,2)
        plot(broadbandSlowWave(IDX==2,1),EMG(IDX==2),'b.')
        hold on
        plot(broadbandSlowWave(EMG>EMGthresh & IDX==1,1),EMG(EMG>EMGthresh & IDX==1),'k.')
        plot(broadbandSlowWave(EMG<EMGthresh & IDX==1|IDX==3,1),EMG(EMG<EMGthresh & IDX==1|IDX==3),'.','Color',0.8*[1 1 1])
        plot(swthresh*[1 1],get(gca,'ylim'),'r','LineWidth',1)
        plot(swthresh*[0 1],EMGthresh*[1 1],'r','LineWidth',1)
        xlabel('Broadband SW');ylabel('EMG')
	subplot(2,2,4)
        %scatter(thratio(SWStimes==0,1),EMG(SWStimes==0,1),3,IDX(SWStimes==0),'filled')
        plot(thratio(NREMtimes==0 & IDX==1,1),EMG(NREMtimes==0 & IDX==1,1),'k.')
        hold on
        plot(thratio(NREMtimes==0 & IDX==3,1),EMG(NREMtimes==0 & IDX==3,1),'r.')
        xlabel('Narrowband Theta');ylabel('EMG')
        plot(THthresh*[1 1],EMGthresh*[0 1],'r','LineWidth',1)
        plot([0 1],EMGthresh*[1 1],'r','LineWidth',1)

saveas(gcf,[figloc,recordingname,'_SSCluster2D'],'jpeg')
%saveas(gcf,['/Users/dlevenstein/Code Library/SleepScoreDevelopment/StateScoreFigures/','ThetaEMGExample'],'jpeg')
%% Figure: Clustering
colormat = [[0 0 0];[0 0 1];[1 0 0]];
coloridx = colormat(IDX,:);

figure
    subplot(1,3,[2,3])
        hold all
        scatter3(broadbandSlowWave,thratio,EMG,2,coloridx,'filled')
        %rotate3d
        view(133.7,18.8);
        grid on
        xlabel('Broadband SW');ylabel('Narrowband Theta');zlabel('EMG')
      
	subplot(3,3,1)
        hold on
        bar(swhistbins,swhist,'FaceColor','none','barwidth',0.9,'linewidth',2)
        plot([swthresh swthresh],[0 max(swhist)],'r')
        xlabel('Slow Wave')
        title('Step 1: Broadband for NREM')
	subplot(3,3,4)
        hold on
        bar(EMGhistbins,EMGhist,'FaceColor','none','barwidth',0.9,'linewidth',2)
        plot([EMGthresh EMGthresh],[0 max(EMGhist)],'r')
        xlabel('EMG')
        title('Step 2: EMG for Muscle Tone')
	subplot(3,3,7)
        hold on
        bar(THhistbins,THhist,'FaceColor','none','barwidth',0.9,'linewidth',2)
        plot([THthresh THthresh],[0 max(THhist)],'r')
        xlabel('Theta')
        title('Step 3: Theta for REM')
        
	saveas(gcf,[figloc,recordingname,'_SSCluster3D'],'jpeg')
%saveas(gcf,['/Users/dlevenstein/Code Library/SleepScoreDevelopment/StateScoreFigures/','clust'],'jpeg')    
  %% Figure: Duration Distributions
%   Wstateintervalss = stateintervals{1};
%   Wlengths = Wstateintervalss(:,2)-Wstateintervalss(:,1);
%   Sstateintervalss = stateintervals{2};
%   Slengths = Sstateintervalss(:,2)-Sstateintervalss(:,1);
%   Rstateintervalss = stateintervals{3};
%   Rlengths = Rstateintervalss(:,2)-Rstateintervalss(:,1);
%   
%   figure
%     subplot(2,3,1)
%         hist(log10(Wlengths),10)
%         set(gca,'XTick',0:3)
%         set(gca,'XTickLabel',10.^[0:3])
%         xlabel('Duration (s)')
%         title('Wake stateintervalserval Durations')
%     subplot(2,3,2)
%         hist(log10(Slengths),10)
%         set(gca,'XTick',0:3)
%         set(gca,'XTickLabel',10.^[0:3])
%         xlabel('Duration (s)')
%         title('SWS stateintervalserval Durations')
%     subplot(2,3,3)
%         hist(log10(Rlengths),10)
%         set(gca,'XTick',0:3)
%         set(gca,'XTickLabel',10.^[0:3])
%         xlabel('Duration (s)')
%         title('REM stateintervalserval Durations')
%     subplot(2,3,4)
%         plot(log10(Wlengths(1:end-1)),log10(Wlengths(2:end)),'.')
%         set(gca,'YTick',0:3)
%         set(gca,'YTickLabel',10.^[0:3])
%         set(gca,'XTick',0:3)
%         set(gca,'XTickLabel',10.^[0:3])
%         xlabel('stateintervalserval n Duration')
%         ylabel('stateintervalserval n+1 Duration')
%         title('Wake stateintervalserval Durations')
%     subplot(2,3,5)
%         plot(log10(Slengths(1:end-1)),log10(Slengths(2:end)),'.')
%         set(gca,'YTick',0:3)
%         set(gca,'YTickLabel',10.^[0:3])
%         set(gca,'XTick',0:3)
%         set(gca,'XTickLabel',10.^[0:3])
%         xlabel('stateintervalserval n Duration')
%         ylabel('stateintervalserval n+1 Duration')
%         title('SWS stateintervalserval Durations')
%     subplot(2,3,6)
%         plot(log10(Rlengths(1:end-1)),log10(Rlengths(2:end)),'.')
%         set(gca,'YTick',0:3)
%         set(gca,'YTickLabel',10.^[0:3])
%         set(gca,'XTick',0:3)
%         set(gca,'XTickLabel',10.^[0:3])
%         xlabel('stateintervalserval n Duration')
%         ylabel('stateintervalserval n+1 Duration')
%         title('REM stateintervalserval Durations')
%         
%         saveas(gcf,[figloc,recordingname,'_stateintervalsdur'],'jpeg')
% %         
%     

 end