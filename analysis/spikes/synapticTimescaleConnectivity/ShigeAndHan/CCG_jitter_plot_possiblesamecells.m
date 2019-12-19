function CCG_jitter_plot_possiblesamecells(ccgjitteroutput)
% Takes output of CCG_jitter_all and plots the cells that function found
% might actually be the same cell.  It plots their CCGs in figures of 5x5
% aaxes each and also labels which shanks they were on and what their ID is
% on that shank, so user may later review.
% Brendon Watson 2012

% 
% %% Retrieve some basic ID info about the cells that seem overlapping
% psc = ccgjitteroutput(1).PossibleSameCells;
% shanks = ccgjitteroutput(1).iEleClu(:,2);
% clus = ccgjitteroutput(1).iEleClu(:,3);
% psc_ids = [shanks(psc(:,1)),clus(psc(:,1)),shanks(psc(:,2)),clus(psc(:,2))];

%% Go thru each possible same cell, gather data about each cell
for x = 1:size(ccgjitteroutput(1).PossibleSameCells,1)

    a = ccgjitteroutput(1).PossibleSameCells(x,1);
    b = ccgjitteroutput(1).PossibleSameCells(x,2);
    
    tR = ccgjitteroutput(a,b).tR;
    ccgR = ccgjitteroutput(a,b).ccgR(:,1,2);
    ccgjm = mean(ccgjitteroutput(a,b).ccgjMtx,2);
    ccgjptMax = ccgjitteroutput(a,b).ccgjstats.pointwiseMax;
    ccgjptMin = ccgjitteroutput(a,b).ccgjstats.pointwiseMin;
    ccgjgbMax = ccgjitteroutput(a,b).ccgjstats.globalMax;
    ccgjgbMin = ccgjitteroutput(a,b).ccgjstats.globalMin;
    GSPExc = ccgjitteroutput(a,b).GSPExc;
    GSPInh = ccgjitteroutput(a,b).GSPInh;
 
%% setting up plotting
    h = figure(ceil(x/25))
    figname(['CCGPossibleSameCellFig#',num2str(x)])
    set(h,'position',[800 40 1100 925],'name',figname);
     
    mod25=mod(x,25);
    if mod25==0
        mod25=25;
    end
    subplot(5,5,mod25)
    
%% Plot
    bar(tR,ccgR,'facecolor','k','edgecolor','k')
    line(tR,ccgjm,'linestyle','--','color','b')
    line(tR,ccgjptMax,'linestyle','--','color','r')
    line(tR,ccgjptMin,'linestyle','--','color','r')
    line(tR,ccgjgbMax,'linestyle','--','color','m')
    line(tR,ccgjgbMin,'linestyle','--','color','m')
    set(gca,'XLim',[min(tR),max(tR)])
    hold on

    TitleText = ['Cells ',num2str(a),' and ', num2str(b),'.'];
    
%% Put marker on significantly positive bin
    findE = find(GSPExc); %find any time bins with significance for excitation from prior ccg_jitter.m
    timeE = tR(findE);%find times of those positive bins
    for c = 1:length(timeE);
      plot(tR(findE(c)),ccgR(findE(c)),'marker','*','color','g')
    end
    title(TitleText);

end

