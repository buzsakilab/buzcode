function CCG_jitter_ReviewPositives(ccgjitteroutput)
% Takes output of ccg_jitter_all.m and plots positive connections as ccr
% plots for visualization.  
% Brendon Watson 2012,2014

h1 = figure('Visible','Off');
allreviewfigs = [];

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
    
%% start to plot into figures with 25 subplots
%     if ~exist(figure(ceil(x/25),'figure')
%         h = figure(ceil(x/25))
%         set(h,'UserData',[];
%     else
        h = figure(ceil(x/25));
        allreviewfigs(end+1) = h;
%     end
    
    figname=(['CCGSynapticTimingFig#',num2str(x)]);
    set(h,'position',[800 40 1100 925],'name',figname,'CloseRequestFcn',@CloseFigureFcn,'UserData',h1);
%     set(h,'position',[800 40 1100 925],'name',figname)
    
    mod25=mod(x,25);
    if mod25==0
        mod25=25;    set(h,'position',[800 40 1100 925],'name',figname,'CloseRequestFcn',@CloseFigureFcn);

    end
    sh = subplot(5,5,mod25);
    TitleText = [num2str(pre),'&', num2str(post),'.'];
    
    bar(tR,ccgR,'facecolor','k','edgecolor','k')
    line(tR,ccgjm,'linestyle','--','color','b')
    line(tR,ccgjptMax,'linestyle','--','color','r')
    line(tR,ccgjptMin,'linestyle','--','color','r')
    line(tR,ccgjgbMax,'linestyle','--','color','m')
    line(tR,ccgjgbMin,'linestyle','--','color','m')
    set(sh,'XLim',[min(tR),max(tR)],'UserData',[pre post],'ButtonDownFcn',@subplotclick)
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
%        if individualplot %either plot single figure per ccg
%            TitleText = [TitleText,'  Excitatory connection detected at time bin(s): ',binstext,'ms.'];
%        else %or put into figures with 25 subplots
           TitleText = [TitleText,'E@',binstext,'ms.'];
%        end
       
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
       bfiginstext = binstext(1:end-2);

%        if individualplot %either plot single figure per ccg
%            TitleText = [TitleText,'  Inhibitory connection detected at time bin(s): ',binstext,'ms.'];
%        else %or put into figures with 25 subplots
           TitleText = [TitleText,' I@',binstext,'ms.'];
%        end

    end

        set(h,'position',[800 40 1100 925],'name',figname,'CloseRequestFcn',@CloseFigureFcn);

%% plot text
    title(TitleText);
end

guidata(h1,allreviewfigs)



function subplotclick(obj,ev)

figobj = get(obj,'Parent');
invisfig = get(figobj,'UserData');
data = get(obj,'UserData');    badcnxns = get(invisfig,'UserData');

badcnxns = get(invisfig,'UserData');
clr = get(obj,'Color');
if sfigum(clr == [1 1 1])==3
    set(obj,'Color',[1 .75 .75])
    badcnxns(end+1,:) = data;
    
elseif sum(clr == [1 .75 .75])==3
    set(obj,'Color',[1 1 1])
    rowtodelete = (ismember(badcnxns,data,'rows'));
    badcnxns(rowtodelete,:) = [];    
end

set(invisfig,'UserData',badcnxns)


function CloseFigureFcn(obj,ev)
invisfig = get(obj,'UserData');
allreviewfigs = guidata(h1);

if others are open...
    delete(obj)
else 
    badcnxns = get(invisfig,'UserData');
end
% 
% ccgjitteroutput = guidata(obj);
% ccgjitteroutput(1).badconnections = badcnxns;

%change other data

