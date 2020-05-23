function    [SleepState] = bz_ThetaStates(basePath,makePlot)
% [SleepState] = bz_ThetaStates(basePath)

% Takes the output of SleepScoreMaster and separtes WAKEstate into
% WAKEtheta and WAKEnontheta by applying the same threshold calculated form
% REM sleep. This is not an ideal way and needs to be improved. 

% F.Sharif, M Valero, AntonioFR still working on it... 5/2020  

if nargin < 2
   makePlot = 1;
end

%% %%%%%%%%%%%
            [~,recordingname] = fileparts(basePath);
            savefolder =basePath;
            bz_sleepstatepath = fullfile(savefolder,[recordingname,'.SleepState.states.mat']);
            load([recordingname '.SleepState.states.mat'])
 
            thratio = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.thratio;
            ththr = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.histsandthreshs.THthresh;
            emg = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.EMG;
            emgthr = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.histsandthreshs.EMGthresh;
            ts = SleepState.detectorinfo.detectionparms.SleepScoreMetrics.t_clus;
            theta_run = thratio>ththr & emg>emgthr;
            Theta_NDX=find(thratio>ththr & emg>emgthr);
            
            A(:,1)=SleepState.idx.states;
            A(Theta_NDX,2)=7;  
            Non_thetaNDX=find(A(:,1)==1& A(:,2)==0);       
            A(Non_thetaNDX,2)=9;  
            
            SleepState.idx.theta_states.states=A(:,2);
            SleepState.idx.theta_states.timestamps=SleepState.idx.timestamps;
            SleepState.idx.theta_states.statenames{7} = 'Theta';
            SleepState.idx.theta_states.statenames{9} = 'Non_theta';
            [INT] = bz_IDXtoINT(SleepState.idx.theta_states);
            SleepState.ints.WAKEtheta=INT.Thetastate;
            SleepState.ints.WAKEnontheta=INT.Non_thetastate;
            save(bz_sleepstatepath,'SleepState');  
            
%%
if makePlot
    
    load([recordingname '.sessionInfo.mat'])

    figure ('position',[0 700 2000 200])
    alpha=.05; Beta=0.09;a=.07; b=.14;
    
    ya=1;
    M=SleepState.ints.WAKEstate;
    if isempty(M)==0
    plot([M(:,1) M(:,2)]/3600,[ya ya],'linewidth',4,'color',[11 102 35]/256)
    dim = [alpha Beta a b];
    str = '--Wake';Cr=[11 102 35]/256;
    annotation('textbox',dim,'String',str,'color',Cr,'FontSize',12,'EdgeColor',Cr,'linewidth',3,'BackgroundColor',Cr,'FaceAlpha',0.2);
    hold on
    end

    ya=2;
    M=SleepState.ints.WAKEtheta;
    if isempty(M)==0
    plot([M(:,1) M(:,2)]/3600,[ya ya],'linewidth',4,'color',[76 187 23]/256)
     dim = [alpha Beta+.15 a b];
    str = '--Thetastate';Cr=[76 187 23]/256;
    annotation('textbox',dim,'String',str,'color',Cr,'FontSize',12,'EdgeColor',Cr,'linewidth',3,'BackgroundColor',Cr,'FaceAlpha',0.2);
    hold on
    end

    ya=3;
    M=SleepState.ints.WAKEnontheta;
    if isempty(M)==0
    plot([M(:,1) M(:,2)]/3600,[ya ya],'linewidth',4,'color',[112 130 56]/256)
    dim = [alpha Beta+.3 a b];
    str = '--Non-thetastate';Cr=[112 130 56]/256;
    annotation('textbox',dim,'String',str,'color',Cr,'FontSize',12,'EdgeColor',Cr,'linewidth',3,'BackgroundColor',Cr,'FaceAlpha',0.2);
    hold on
    end

    ya=6;
    M=SleepState.ints.NREMstate;
    if isempty(M)==0
    plot([M(:,1) M(:,2)]/3600,[ya ya],'linewidth',4,'color','r')
    dim = [alpha Beta+.45 a b];
    str = '--NREMstate';Cr='r';
    annotation('textbox',dim,'String',str,'color',Cr,'FontSize',12,'EdgeColor',Cr,'linewidth',3,'BackgroundColor',Cr,'FaceAlpha',0.2);
    hold on
    end

    ya=7;
    M=SleepState.ints.REMstate;
    if isempty(M)==0
    plot([M(:,1) M(:,2)]/3600,[ya ya],'linewidth',4,'color','b')
    dim = [alpha Beta+.6 a b];
    str = '--REMstate';Cr='b';
    annotation('textbox',dim,'String',str,'color',Cr,'FontSize',12,'EdgeColor',Cr,'linewidth',3,'BackgroundColor',Cr,'FaceAlpha',0.2);
    hold on
    end


    set(gcf,'Color','w')
    xlabel('Recording time (hr)','FontSize',12)
    title([recordingname ' - Recording States'],'FontSize',12)
    yticks([])
    box off    
    
    saveas(gcf,strcat(savefolder,filesep,'StateScoreFigures',filesep,recordingname,'_thetaStates'),'jpg');
    
end

