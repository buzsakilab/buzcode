function [ INT, IDX, t_IDX,PC1weights,PC1expvar,broadbandSlowWave,thratio,EMG,t_FFT ,badtimes, reclength] = ClusterStates(swLFP,thLFP,EMG,sf_LFP,sf_EMG,figloc,recordingname)
%StateID(LFP,thLFP,EMG,sf_LFP,sf_EMG,figloc,WSEpisodes)
%   Detailed explanation goes here
%
%
%
%Dependencies: IDXtoINT_ss, bz_INTtoIDX
%
%Last Updated: 1/31/16
%DLevenstein




%% Min Win Parameters (s)
%if exist(
minSWS = 6;
minWnexttoREM = 6;
minWinREM = 6;       
minREMinW = 6;
minREM = 6;
minWAKE = 6;


%% Downsample and filter
%Make Downsample to niquest frequency

if sf_LFP == 1250
    downsamplefactor = 5;
elseif sf_LFP == 250
    downsamplefactor = 1;
elseif sf_LFP == 1000
    downsamplefactor = 4;
else
    display('sf not recognized... if only you made this able to set its own downsample...')
end
swLFP = downsample(swLFP,downsamplefactor);
thLFP = downsample(thLFP,downsamplefactor);
sf_LFP = sf_LFP/downsamplefactor;


%filtbounds = [0.5 120];
%display(['Filtering ',num2str(filtbounds(1)),'-',num2str(filtbounds(2)),' Hz...']);
%LFP = FiltNPhase(LFP, filtbounds, sf_LFP );


%% Calculate Spectrogram
%display('FFT Spectrum for Broadband LFP')

freqlist = logspace(0,2,100);
%freqlist = linspace(0.5,55.5,111);
window = 10;
noverlap = 9;
window = window*sf_LFP;
noverlap = noverlap*sf_LFP;
[FFTspec,FFTfreqs,t_FFT] = spectrogram(swLFP,window,noverlap,freqlist,sf_LFP);
FFTspec = abs(FFTspec);
[zFFTspec,mu,sig] = zscore(log10(FFTspec)');

    %% Remove transients before calculating SW histogram
    %this should be it's own whole section - removing/detecting transients
totz = zscore(abs(sum(zFFTspec')));
badtimes = find(totz>5);
zFFTspec(badtimes,:) = 0;

%% PCA for Broadband Slow Wave
 [COEFF, SCORE, ~, ~, EXPLAINED] = pca(zFFTspec);
 % broadbandSlowWave = SCORE(:,1);
PC1weights = COEFF(:,1);
PC1expvar = EXPLAINED(1);
 
%% Set Broadband filter weights for Slow Wave
load('SWweights.mat')
assert(isequal(freqlist,SWfreqlist), 'spectrogram freqs.  are not what they should be...')
broadbandSlowWave = zFFTspec*SWweights';
 
%% Smooth and 0-1 normalize
smoothfact = 10; %units of si_FFT
thsmoothfact = 10; %used to be 15

broadbandSlowWave = smooth(broadbandSlowWave,smoothfact);
broadbandSlowWave = (broadbandSlowWave-min(broadbandSlowWave))./max(broadbandSlowWave-min(broadbandSlowWave));


 
%% Calculate theta
%display('FFT Spectrum for Theta')

% %NarrowbandTheta
f_all = [2 20];
f_theta = [5 10];
freqlist = logspace(log10(f_all(1)),log10(f_all(2)),100);


[thFFTspec,thFFTfreqs] = spectrogram(thLFP,window,noverlap,freqlist,sf_LFP);
thFFTspec = (abs(thFFTspec));
[~,mu_th,sig_th] = zscore(log10(thFFTspec)');

thfreqs = find(thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
allpower = sum((thFFTspec),1);
thpower = sum((thFFTspec(thfreqs,:)),1);

thratio = thpower./allpower;    %Narrowband Theta
thratio = smooth(thratio,thsmoothfact);
thratio = (thratio-min(thratio))./max(thratio-min(thratio));
 
%% EMG
dtEMG = 1/sf_EMG;
t_EMG = (1:length(EMG))*dtEMG;
EMG = smooth(EMG,smoothfact/dtEMG);
EMG = (EMG-min(EMG))./max(EMG-min(EMG));

reclength = round(t_EMG(end));

%downsample to FFT time points;
[~,t_intersect] = intersect(t_EMG,t_FFT);
EMG = EMG(t_intersect);
t_EMG = t_EMG(t_intersect);



%% Divide PC1 for SWS
numpeaks = 1;
numbins = 12;
%numbins = 12; %for Poster...
while numpeaks ~=2
    [pcahist,histbins]= hist(broadbandSlowWave,numbins);
    
    [PKS,LOCS] = findpeaks_SleepScore(pcahist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end


betweenpeaks = histbins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks_SleepScore(-pcahist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

thresh = betweenpeaks(diploc);

%Set transients to wake state
broadbandSlowWave(badtimes,1)=histbins(LOCS(1));
 
 
%SWS time points
NREMtimes = (broadbandSlowWave >thresh);


%% Then Divide EMG
numpeaks = 1;
numbins = 12;
while numpeaks ~=2
    [EMGhist,EMGhistbins]= hist(EMG(NREMtimes==0),numbins);
    %[EMGhist,EMGhistbins]= hist(EMG,numbins);

    [PKS,LOCS] = findpeaks_SleepScore([0 EMGhist],'NPeaks',2);
    LOCS = sort(LOCS)-1;
    numbins = numbins+1;
    numpeaks = length(LOCS);
    
    if numpeaks ==100
        display('Something is wrong with your EMG')
        return
    end
end

betweenpeaks = EMGhistbins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks_SleepScore(-EMGhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

EMGthresh = betweenpeaks(diploc);

MOVtimes = (broadbandSlowWave<thresh & EMG>EMGthresh);


%% Then Divide Theta
numpeaks = 1;
numbins = 12;
while numpeaks ~=2 && numbins <=25
    %[THhist,THhistbins]= hist(thratio(SWStimes==0 & MOVtimes==0),numbins);
    [THhist,THhistbins]= hist(thratio(MOVtimes==0),numbins);

    [PKS,LOCS] = findpeaks_SleepScore(THhist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

numbins = 12;
%numbins = 15; %for Poster...
while numpeaks ~=2 && numbins <=25
    [THhist,THhistbins]= hist(thratio(NREMtimes==0 & MOVtimes==0),numbins);

    [PKS,LOCS] = findpeaks_SleepScore(THhist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

if length(PKS)==2
    betweenpeaks = THhistbins(LOCS(1):LOCS(2));
    [dip,diploc] = findpeaks_SleepScore(-THhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

    THthresh = betweenpeaks(diploc);

    REMtimes = (broadbandSlowWave<thresh & EMG<EMGthresh & thratio>THthresh);
else
    THthresh = 0;
    REMtimes =(broadbandSlowWave<thresh & EMG<EMGthresh);
end

%%
%OLD:
%Index Vector: SWS=2, REM=3, MOV=6, NonMOV=1.   
%(Separate MOV for REM, then join later)
%IDX = SWStimes+2*REMtimes+5*MOVtimes+1;

%NEW: No separation of MOV and NonMOV WAKE
%Index Vector: NREM=2, REM=3, WAKE=1. 
IDX = NREMtimes+2*REMtimes+1;

%Start/end offset due to FFT


%% Minimum Interuptions
INT = IDXtoINT_ss(IDX,3);


%Make the following repeated chunks of code into a single function.

%SWS  (to NonMOV)
Sints = INT{2};
Slengths = Sints(:,2)-Sints(:,1);
shortSints = {Sints(find(Slengths<=minSWS),:)};
shortSidx = bz_INTtoIDX(shortSints,'length',length(IDX));
%Change Short SWS to Wake
IDX(shortSidx==1) = 1;   
INT = IDXtoINT_ss(IDX,3);

%NonMOV next to REM   (to REM)
Wints = INT{1};
trans = (diff(IDX)); %All State Transitions
WRtrans = find((trans)==-2)+1;  %Just transitions between WAKE, REM
%Convert to interval indices
[~,WRtransON] = intersect(Wints(:,1),WRtrans);
WRtrans = find((trans)==2);
[~,WRtransOFF] = intersect(Wints(:,2),WRtrans);
WRtrans = union(WRtransON,WRtransOFF); %On or offset are RW
%Find WAKE intervals that border REM and are less than min
Wlengths = Wints(:,2)-Wints(:,1);
shortWRints = find(Wlengths(WRtrans)<=minWnexttoREM);
shortWRints = WRtrans(shortWRints);
shortWRints = {Wints(shortWRints,:)};
shortWRidx = bz_INTtoIDX(shortWRints,'length',length(IDX));
%Convert wake to rem
IDX(shortWRidx==1) = 3;
INT = IDXtoINT_ss(IDX,3);


%NonMOV in REM   (to REM)
Wints = INT{1};
trans = (diff(IDX)); %All State Transitions
WRtrans = find((trans)==-2)+1;  %Just transitions between WAKE, REM
%Convert to interval indices
[~,WRtransON] = intersect(Wints(:,1),WRtrans);
WRtrans = find((trans)==2);
[~,WRtransOFF] = intersect(Wints(:,2),WRtrans);
WRtrans = intersect(WRtransON,WRtransOFF); %Both onset and offset are RW
%Find WAKE intervals that border REM and are less than min
Wlengths = Wints(:,2)-Wints(:,1);
shortWRints = find(Wlengths(WRtrans)<=minWinREM);
shortWRints = WRtrans(shortWRints);
shortWRints = {Wints(shortWRints,:)};
shortWRidx = bz_INTtoIDX(shortWRints,'length',length(IDX));
%Convert wake to rem
IDX(shortWRidx==1) = 3;
IDX(IDX==6) = 1; %Convert NonMOV to WAKE
INT = IDXtoINT_ss(IDX,3);


%REM in WAKE   (to WAKE)
Rints = INT{3};
trans = (diff(IDX)); %All State Transitions
WRtrans = find((trans)==2)+1;  %Just transitions between WAKE, REM
%Convert to interval indices
[~,WRtransON] = intersect(Rints(:,1),WRtrans);
WRtrans = find((trans)==-2);
[~,WRtransOFF] = intersect(Rints(:,2),WRtrans);
WRtrans = intersect(WRtransON,WRtransOFF); %Both onset and offset are RW
%Find WAKE intervals that border REM and are less than min
Rlengths = Rints(:,2)-Rints(:,1);
shortWRints = find(Rlengths(WRtrans)<=minREMinW);
shortWRints = WRtrans(shortWRints);
shortWRints = {Rints(shortWRints,:)};
shortWRidx = bz_INTtoIDX(shortWRints,'length',length(IDX));
%Convert REM to WAKE
IDX(shortWRidx==1) = 1;
INT = IDXtoINT_ss(IDX,3);

%REM (only applies to REM in the middle of SWS)    (to WAKE)
Rints = INT{3};
Rlengths = Rints(:,2)-Rints(:,1);
shortRints = {Rints(find(Rlengths<=minREM),:)};
shortRidx = bz_INTtoIDX(shortRints,'length',length(IDX));

IDX(shortRidx==1) = 1;
INT = IDXtoINT_ss(IDX,3);


%WAKE   (to SWS)     essentiall a minimum MA time
Wints = INT{1};
Wlengths = Wints(:,2)-Wints(:,1);
shortWints = {Wints(find(Wlengths<=minWAKE),:)};
shortWidx = bz_INTtoIDX(shortWints,'length',length(IDX));
IDX(shortWidx==1) = 2;

INT = IDXtoINT_ss(IDX,3);

%SWS  (to NonMOV)
Sints = INT{2};
Slengths = Sints(:,2)-Sints(:,1);
shortSints = {Sints(find(Slengths<=minSWS),:)};
shortSidx = bz_INTtoIDX(shortSints,'length',length(IDX));
%Change Short SWS to Wake
IDX(shortSidx==1) = 1;   
INT = IDXtoINT_ss(IDX,3);




%% Pad time to match recording time
offset = t_FFT(1)-1;

INT = cellfun(@(x) x+offset,INT,'UniformOutput',false);


 %% Figure
 
 if exist('figloc','var')
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
        

        
%% 
figure

    subplot(2,3,1)
        scatter(broadbandSlowWave,thratio,3,IDX,'filled')
        xlabel('Broadband SW');ylabel('Narrowband Theta')
    subplot(2,3,2)
        scatter(broadbandSlowWave,EMG,3,IDX,'filled')
        xlabel('Broadband SW');ylabel('EMG')
    subplot(2,3,3)
        scatter(thratio,EMG,3,IDX,'filled')
        xlabel('Narrowband Theta');ylabel('EMG')

    subplot(2,3,4)
        scatter(broadbandSlowWave(NREMtimes==0,1),thratio(NREMtimes==0),3,IDX(NREMtimes==0),'filled')
        xlabel('Broadband SW');ylabel('Narrowband Theta')
    subplot(2,3,5)
        scatter(broadbandSlowWave(NREMtimes==0,1),EMG(NREMtimes==0),3,IDX(NREMtimes==0),'filled')
        xlabel('Broadband SW');ylabel('EMG')
        title('non-nonREM only')
    subplot(2,3,6)
        %scatter(thratio(SWStimes==0,1),EMG(SWStimes==0,1),3,IDX(SWStimes==0),'filled')
        plot(thratio(NREMtimes==0 & IDX==1,1),EMG(NREMtimes==0 & IDX==1,1),'k.')
        hold on
        plot(thratio(NREMtimes==0 & IDX==3,1),EMG(NREMtimes==0 & IDX==3,1),'r.')
        xlabel('Narrowband Theta');ylabel('EMG')

%% Figure: Split REM/Arousal  
figure
	subplot(3,2,1)
        hold on
        bar(histbins(histbins>thresh),pcahist(histbins>thresh),'FaceColor','b','barwidth',0.9,'linewidth',1)
        bar(histbins(histbins<=thresh),pcahist(histbins<=thresh),'FaceColor',0.9*[1 1 1],'barwidth',0.9,'linewidth',1)
        plot([thresh thresh],[0 max(pcahist)],'r','LineWidth',1)
        xlabel('PC 1')
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
        plot(thresh*[1 1],get(gca,'ylim'),'r','LineWidth',1)
        plot(thresh*[0 1],EMGthresh*[1 1],'r','LineWidth',1)
        xlabel('Broadband SW');ylabel('EMG')
	subplot(2,2,4)
        %scatter(thratio(SWStimes==0,1),EMG(SWStimes==0,1),3,IDX(SWStimes==0),'filled')
        plot(thratio(NREMtimes==0 & IDX==1,1),EMG(NREMtimes==0 & IDX==1,1),'k.')
        hold on
        plot(thratio(NREMtimes==0 & IDX==3,1),EMG(NREMtimes==0 & IDX==3,1),'r.')
        xlabel('Narrowband Theta');ylabel('EMG')
        plot(THthresh*[1 1],EMGthresh*[0 1],'r','LineWidth',1)
        plot([0 1],EMGthresh*[1 1],'r','LineWidth',1)

saveas(gcf,[figloc,recordingname,'_clust2'],'jpeg')
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
        bar(histbins,pcahist,'FaceColor','none','barwidth',0.9,'linewidth',2)
        plot([thresh thresh],[0 max(pcahist)],'r')
        xlabel('PC 1')
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
        
	saveas(gcf,[figloc,recordingname,'_clust'],'jpeg')
%saveas(gcf,['/Users/dlevenstein/Code Library/SleepScoreDevelopment/StateScoreFigures/','clust'],'jpeg')    
  %% Figure: Duration Distributions
%   Wints = INT{1};
%   Wlengths = Wints(:,2)-Wints(:,1);
%   Sints = INT{2};
%   Slengths = Sints(:,2)-Sints(:,1);
%   Rints = INT{3};
%   Rlengths = Rints(:,2)-Rints(:,1);
%   
%   figure
%     subplot(2,3,1)
%         hist(log10(Wlengths),10)
%         set(gca,'XTick',0:3)
%         set(gca,'XTickLabel',10.^[0:3])
%         xlabel('Duration (s)')
%         title('Wake Interval Durations')
%     subplot(2,3,2)
%         hist(log10(Slengths),10)
%         set(gca,'XTick',0:3)
%         set(gca,'XTickLabel',10.^[0:3])
%         xlabel('Duration (s)')
%         title('SWS Interval Durations')
%     subplot(2,3,3)
%         hist(log10(Rlengths),10)
%         set(gca,'XTick',0:3)
%         set(gca,'XTickLabel',10.^[0:3])
%         xlabel('Duration (s)')
%         title('REM Interval Durations')
%     subplot(2,3,4)
%         plot(log10(Wlengths(1:end-1)),log10(Wlengths(2:end)),'.')
%         set(gca,'YTick',0:3)
%         set(gca,'YTickLabel',10.^[0:3])
%         set(gca,'XTick',0:3)
%         set(gca,'XTickLabel',10.^[0:3])
%         xlabel('Interval n Duration')
%         ylabel('Interval n+1 Duration')
%         title('Wake Interval Durations')
%     subplot(2,3,5)
%         plot(log10(Slengths(1:end-1)),log10(Slengths(2:end)),'.')
%         set(gca,'YTick',0:3)
%         set(gca,'YTickLabel',10.^[0:3])
%         set(gca,'XTick',0:3)
%         set(gca,'XTickLabel',10.^[0:3])
%         xlabel('Interval n Duration')
%         ylabel('Interval n+1 Duration')
%         title('SWS Interval Durations')
%     subplot(2,3,6)
%         plot(log10(Rlengths(1:end-1)),log10(Rlengths(2:end)),'.')
%         set(gca,'YTick',0:3)
%         set(gca,'YTickLabel',10.^[0:3])
%         set(gca,'XTick',0:3)
%         set(gca,'XTickLabel',10.^[0:3])
%         xlabel('Interval n Duration')
%         ylabel('Interval n+1 Duration')
%         title('REM Interval Durations')
%         
%         saveas(gcf,[figloc,recordingname,'_intdur'],'jpeg')
%         
    

 end
 
IDX = bz_INTtoIDX(INT,'length',reclength);
t_IDX = 1:length(IDX);
IDX = IDX';

end

