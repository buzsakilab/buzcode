function ClusterStates_MakeFigure(INT,IDX,figloc,FFTfreqs,FFTspec,thFFTfreqs,thFFTspec,t_FFT,recordingname,broadbandSlowWave,thratio,EMG,t_EMG) 

%% Figure
[zFFTspec,mu,sig] = zscore(log10(FFTspec)');
[~,mu_th,sig_th] = zscore(log10(thFFTspec)');

 viewwin  =[t_FFT(1) t_FFT(end)];
 %viewwin  =[32000 34000];
%viewwin=[9000 11000];
figure
	subplot(8,1,[1:2])
        imagesc(t_FFT,log2(FFTfreqs),log10(FFTspec))
        axis xy
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        caxis([3.5 6.5])
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        xlim(viewwin)
        colorbar('east')
        ylim([log2(FFTfreqs(1)) log2(FFTfreqs(end))+0.2])
        set(gca,'XTickLabel',{})
        ylabel({'swLFP','f (Hz)'})
        title([recordingname,': State Scoring Results']);
	subplot(8,1,3)
        imagesc(t_FFT,log2(thFFTfreqs),log10(thFFTspec))
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
        %plot(t_FFT,-IDX,'LineWidth',2)
        hold on
        plot(INT{1}',-1*ones(size(INT{1}))','k','LineWidth',8)
        plot(INT{2}',-2*ones(size(INT{2}))','b','LineWidth',8)
        plot(INT{3}',-3*ones(size(INT{3}))','r','LineWidth',8)
        xlim(viewwin)
        ylim([-4 0])
        set(gca,'YTick',[-3:-1])
        set(gca,'YTickLabel',{'REM','SWS','Wake/MA'})
        set(gca,'XTickLabel',{})
        
   	subplot(6,1,4)
        hold on
        plot(t_FFT,broadbandSlowWave,'k')
        %plot(synchtimes',thresh*ones(size(synchtimes))','r')
        ylabel('SW')
        xlim([t_FFT(1) t_FFT(end)])
        xlim(viewwin)
        set(gca,'XTickLabel',{})
        
   	subplot(6,1,5)
        hold on
        plot(t_FFT,thratio,'k')
        %plot(synchtimes',thresh*ones(size(synchtimes))','r')
        ylabel('Theta')
        xlim([t_FFT(1) t_FFT(end)])
        xlim(viewwin)
        set(gca,'XTickLabel',{})
        
   	subplot(6,1,6)
        hold on
        plot(t_EMG,EMG,'k')
        %plot(synchtimes',thresh*ones(size(synchtimes))','r')
        ylabel('EMG')
        xlim([t_FFT(1) t_FFT(end)])
        xlim(viewwin)
        xlabel('t (s)')
        
	saveas(gcf,[figloc,recordingname,'_ClusterResults'],'jpeg')
        

        
%% didn't get this far yet as far as importing variables etc...
%% ... gotta keep going, get NREMtimes etc
% figure
% 
%     subplot(2,3,1)
%         scatter(broadbandSlowWave,thratio,3,IDX,'filled')
%         xlabel('Broadband SW');ylabel('Narrowband Theta')
%     subplot(2,3,2)
%         scatter(broadbandSlowWave,EMG,3,IDX,'filled')
%         xlabel('Broadband SW');ylabel('EMG')
%     subplot(2,3,3)
%         scatter(thratio,EMG,3,IDX,'filled')
%         xlabel('Narrowband Theta');ylabel('EMG')
% 
%     subplot(2,3,4)
%         scatter(broadbandSlowWave(NREMtimes==0,1),thratio(NREMtimes==0),3,IDX(NREMtimes==0),'filled')
%         xlabel('Broadband SW');ylabel('Narrowband Theta')
%     subplot(2,3,5)
%         scatter(broadbandSlowWave(NREMtimes==0,1),EMG(NREMtimes==0),3,IDX(NREMtimes==0),'filled')
%         xlabel('Broadband SW');ylabel('EMG')
%         title('non-nonREM only')
%     subplot(2,3,6)
%         %scatter(thratio(SWStimes==0,1),EMG(SWStimes==0,1),3,IDX(SWStimes==0),'filled')
%         plot(thratio(NREMtimes==0 & IDX==1,1),EMG(NREMtimes==0 & IDX==1,1),'k.')
%         hold on
%         plot(thratio(NREMtimes==0 & IDX==3,1),EMG(NREMtimes==0 & IDX==3,1),'r.')
%         xlabel('Narrowband Theta');ylabel('EMG')
% 
% %% Figure: Split REM/Arousal  
% figure
% 	subplot(3,2,1)
%         hold on
%         bar(histbins(histbins>swthresh),pcahist(histbins>swthresh),'FaceColor','b','barwidth',0.9,'linewidth',1)
%         bar(histbins(histbins<=swthresh),pcahist(histbins<=swthresh),'FaceColor',0.9*[1 1 1],'barwidth',0.9,'linewidth',1)
%         plot([swthresh swthresh],[0 max(pcahist)],'r','LineWidth',1)
%         xlabel('PC 1')
%         title('Step 1: Broadband for NREM')
%         
% 
% 	subplot(3,2,3)
%         hold on
%         bar(EMGhistbins(EMGhistbins>EMGthresh),EMGhist(EMGhistbins>EMGthresh),'FaceColor','k','barwidth',0.9,'linewidth',1)
%         bar(EMGhistbins(EMGhistbins<=EMGthresh),EMGhist(EMGhistbins<=EMGthresh),'FaceColor',0.9*[1 1 1],'barwidth',0.9,'linewidth',1)
%         plot([EMGthresh EMGthresh],[0 max(EMGhist)],'r','LineWidth',1)
%         xlabel('EMG')
%         title('Step 2: EMG for Muscle Tone')
% 	subplot(3,2,5)
%         hold on
%         bar(THhistbins(THhistbins>=THthresh),THhist(THhistbins>=THthresh),'FaceColor','r','barwidth',0.9,'linewidth',1)
%         bar(THhistbins(THhistbins<THthresh),THhist(THhistbins<THthresh),'FaceColor','k','barwidth',0.9,'linewidth',1)
%         plot([THthresh THthresh],[0 max(THhist)],'r','LineWidth',1)
%         xlabel('Theta')
%         title('Step 3: Theta for REM')
%         
%         
%     subplot(2,2,2)
%         plot(broadbandSlowWave(IDX==2,1),EMG(IDX==2),'b.')
%         hold on
%         plot(broadbandSlowWave(EMG>EMGthresh & IDX==1,1),EMG(EMG>EMGthresh & IDX==1),'k.')
%         plot(broadbandSlowWave(EMG<EMGthresh & IDX==1|IDX==3,1),EMG(EMG<EMGthresh & IDX==1|IDX==3),'.','Color',0.8*[1 1 1])
%         plot(swthresh*[1 1],get(gca,'ylim'),'r','LineWidth',1)
%         plot(swthresh*[0 1],EMGthresh*[1 1],'r','LineWidth',1)
%         xlabel('Broadband SW');ylabel('EMG')
% 	subplot(2,2,4)
%         %scatter(thratio(SWStimes==0,1),EMG(SWStimes==0,1),3,IDX(SWStimes==0),'filled')
%         plot(thratio(NREMtimes==0 & IDX==1,1),EMG(NREMtimes==0 & IDX==1,1),'k.')
%         hold on
%         plot(thratio(NREMtimes==0 & IDX==3,1),EMG(NREMtimes==0 & IDX==3,1),'r.')
%         xlabel('Narrowband Theta');ylabel('EMG')
%         plot(THthresh*[1 1],EMGthresh*[0 1],'r','LineWidth',1)
%         plot([0 1],EMGthresh*[1 1],'r','LineWidth',1)
% 
% saveas(gcf,[figloc,recordingname,'_clust2'],'jpeg')
% %saveas(gcf,['/Users/dlevenstein/Code Library/SleepScoreDevelopment/StateScoreFigures/','ThetaEMGExample'],'jpeg')
% %% Figure: Clustering
% colormat = [[0 0 0];[0 0 1];[1 0 0]];
% coloridx = colormat(IDX,:);
% 
% figure
%     subplot(1,3,[2,3])
%         hold all
%         scatter3(broadbandSlowWave,thratio,EMG,2,coloridx,'filled')
%         %rotate3d
%         view(133.7,18.8);
%         grid on
%         xlabel('Broadband SW');ylabel('Narrowband Theta');zlabel('EMG')
%       
% 	subplot(3,3,1)
%         hold on
%         bar(histbins,pcahist,'FaceColor','none','barwidth',0.9,'linewidth',2)
%         plot([swthresh swthresh],[0 max(pcahist)],'r')
%         xlabel('PC 1')
%         title('Step 1: Broadband for NREM')
% 	subplot(3,3,4)
%         hold on
%         bar(EMGhistbins,EMGhist,'FaceColor','none','barwidth',0.9,'linewidth',2)
%         plot([EMGthresh EMGthresh],[0 max(EMGhist)],'r')
%         xlabel('EMG')
%         title('Step 2: EMG for Muscle Tone')
% 	subplot(3,3,7)
%         hold on
%         bar(THhistbins,THhist,'FaceColor','none','barwidth',0.9,'linewidth',2)
%         plot([THthresh THthresh],[0 max(THhist)],'r')
%         xlabel('Theta')
%         title('Step 3: Theta for REM')
%         
% 	saveas(gcf,[figloc,recordingname,'_clust'],'jpeg')
% %saveas(gcf,['/Users/dlevenstein/Code Library/SleepScoreDevelopment/StateScoreFigures/','clust'],'jpeg')    
%   %% Figure: Duration Distributions
% %   Wints = INT{1};
% %   Wlengths = Wints(:,2)-Wints(:,1);
% %   Sints = INT{2};
% %   Slengths = Sints(:,2)-Sints(:,1);
% %   Rints = INT{3};
% %   Rlengths = Rints(:,2)-Rints(:,1);
% %   
% %   figure
% %     subplot(2,3,1)
% %         hist(log10(Wlengths),10)
% %         set(gca,'XTick',0:3)
% %         set(gca,'XTickLabel',10.^[0:3])
% %         xlabel('Duration (s)')
% %         title('Wake Interval Durations')
% %     subplot(2,3,2)
% %         hist(log10(Slengths),10)
% %         set(gca,'XTick',0:3)
% %         set(gca,'XTickLabel',10.^[0:3])
% %         xlabel('Duration (s)')
% %         title('SWS Interval Durations')
% %     subplot(2,3,3)
% %         hist(log10(Rlengths),10)
% %         set(gca,'XTick',0:3)
% %         set(gca,'XTickLabel',10.^[0:3])
% %         xlabel('Duration (s)')
% %         title('REM Interval Durations')
% %     subplot(2,3,4)
% %         plot(log10(Wlengths(1:end-1)),log10(Wlengths(2:end)),'.')
% %         set(gca,'YTick',0:3)
% %         set(gca,'YTickLabel',10.^[0:3])
% %         set(gca,'XTick',0:3)
% %         set(gca,'XTickLabel',10.^[0:3])
% %         xlabel('Interval n Duration')
% %         ylabel('Interval n+1 Duration')
% %         title('Wake Interval Durations')
% %     subplot(2,3,5)
% %         plot(log10(Slengths(1:end-1)),log10(Slengths(2:end)),'.')
% %         set(gca,'YTick',0:3)
% %         set(gca,'YTickLabel',10.^[0:3])
% %         set(gca,'XTick',0:3)
% %         set(gca,'XTickLabel',10.^[0:3])
% %         xlabel('Interval n Duration')
% %         ylabel('Interval n+1 Duration')
% %         title('SWS Interval Durations')
% %     subplot(2,3,6)
% %         plot(log10(Rlengths(1:end-1)),log10(Rlengths(2:end)),'.')
% %         set(gca,'YTick',0:3)
% %         set(gca,'YTickLabel',10.^[0:3])
% %         set(gca,'XTick',0:3)
% %         set(gca,'XTickLabel',10.^[0:3])
% %         xlabel('Interval n Duration')
% %         ylabel('Interval n+1 Duration')
% %         title('REM Interval Durations')
% %         
% %         saveas(gcf,[figloc,recordingname,'_intdur'],'jpeg')
% %         
%     

 end