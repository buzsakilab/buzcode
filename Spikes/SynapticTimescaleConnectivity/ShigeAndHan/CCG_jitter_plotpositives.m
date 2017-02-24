function CCG_jitter_plotpositives(ccgjitteroutput,individualplot)
% Takes output of ccg_jitter_all.m and plots positive connections as ccr
% plots for visualization.  
% Brendon Watson 2012

if ~exist('individualplot','var');
    individualplot = 0;
end


%% connections are pre-recorded in ccgjitteroutput, gather them up
AllConnections = cat(1,ccgjitteroutput(1).ConnectionsE, ccgjitteroutput(1).ConnectionsI);
ConnectionTypes = cat(1,ones(size(ccgjitteroutput(1).ConnectionsE,1),1),zeros(size(ccgjitteroutput(1).ConnectionsI,1),1));
     % vector of 1's where excitatory, 0's where inhibitory

%% Go thru each connection, gather data for relevant cells
for x = 1:size(AllConnections,1)
    pre = AllConnections(x,1);
    post = AllConnections(x,2);
        
    a = min([pre post]);
    b = max([pre post]);
    
    tR = ccgjitteroutput(a,b).tR;
    ccgR = ccgjitteroutput(a,b).ccgR(:,1,2);
    ccgjm = mean(ccgjitteroutput(a,b).ccgjMtx,2);
    ccgjptMax = ccgjitteroutput(a,b).ccgjstats.pointwiseMax;
    ccgjptMin = ccgjitteroutput(a,b).ccgjstats.pointwiseMin;
    ccgjgbMax = ccgjitteroutput(a,b).ccgjstats.globalMax;
    ccgjgbMin = ccgjitteroutput(a,b).ccgjstats.globalMin;
    GSPExc = ccgjitteroutput(a,b).GSPExc;
    GSPInh = ccgjitteroutput(a,b).GSPInh;
    
%% if a "backwards" connection with 2nd listed cell as presynaptic, then
%% reverse times on variables
    if b==pre
        ccgR = flipud(ccgR);
        ccgjm = flipud(ccgjm);
        ccgjptMax = fliplr(ccgjptMax);
        ccgjptMin = fliplr(ccgjptMin);
        ccgjgbMax = fliplr(ccgjgbMax);
        ccgjgbMin = fliplr(ccgjgbMin);
        GSPExc = fliplr(GSPExc);
        GSPInh = fliplr(GSPInh);
    end
    
%% start to plot
    if individualplot %either plot single figure per ccg
        h = figure;
        TitleText = ['CCG comparisons for cells ',num2str(pre),' and ', num2str(post),'.'];
        set(h,'name',['CCGSynapticTimingCells',num2str(pre),' and ', num2str(post)])
    else %or put into figures with 25 subplots
        h = figure(ceil(x/25));
        figname=(['CCGSynapticTimingFig#',num2str(x)]);
        set(h,'position',[800 40 1100 925],'name',figname);
    
        mod25=mod(x,25);
        if mod25==0
            mod25=25;
        end
        subplot(5,5,mod25)
        TitleText = [num2str(pre),'&', num2str(post),'.'];
    end    
    
    bar(tR,ccgR,'facecolor','k','edgecolor','k')
    line(tR,ccgjm,'linestyle','--','color','b')
    line(tR,ccgjptMax,'linestyle','--','color','r')
    line(tR,ccgjptMin,'linestyle','--','color','r')
    line(tR,ccgjgbMax,'linestyle','--','color','m')
    line(tR,ccgjgbMin,'linestyle','--','color','m')
    set(gca,'XLim',[min(tR),max(tR)])
    hold on
        
%% for excitatory cells, plot significant peaks, add descriptive title
    if ConnectionTypes(x)%if excitatory
       binstext = [];
       findE = find(GSPExc); %find any time bins with significance for excitation from prior ccg_jitter.m
       timeE = tR(findE);%find times of those positive bins
       for c = 1:length(timeE);
          plot(tR(findE(c)),ccgR(findE(c)),'marker','*','color','g')
          binstext = [binstext,num2str(tR(findE(c))),','];
       end
       binstext = binstext(1:end-1);
       if individualplot %either plot single figure per ccg
           TitleText = [TitleText,'  Excitatory connection detected at time bin(s): ',binstext,'ms.'];
       else %or put into figures with 25 subplots
           TitleText = [TitleText,'E@',binstext,'ms.'];
       end
       
    end


%% same for inhibitory
    if ~ConnectionTypes(x)%if inhibitory
       binstext = [];
       findI = find(GSPInh); %find any time bins with significance for excitation from prior ccg_jitter.m
       timeI = tR(findI);%find times of those positive bins
       for c = 1:length(timeI);
          plot(tR(findI(c)),ccgR(findI(c)),'marker','*','color','g')
          binstext = [binstext,num2str(tR(findI(c))),', '];
       end
       binstext = binstext(1:end-2);

       if individualplot %either plot single figure per ccg
           TitleText = [TitleText,'  Inhibitory connection detected at time bin(s): ',binstext,'ms.'];
       else %or put into figures with 25 subplots
           TitleText = [TitleText,' I@',binstext,'ms.'];
       end

    end

    
%% plot text
    title(TitleText);
end