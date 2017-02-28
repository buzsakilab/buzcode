function invisfig = FindSynapse_ReviewOutput(funcsynapses,rawS,outputvarname,allpresyn,badcnxns)
% Takes output of FindSynapse.m and plots positive connections as ccr
% plots for review and clickable deletions of bad connections.
% INPUTS
% funcsynapses - struct array containing population connectivity information,
%            created by Make_FindSynapse_bw and/or 
% rawS - a tsd array of spike times, NOT burst filtered
% outputvarname - the name of the output variable to create in the caller
%            workspace
% allpresyn (optional) - a list of presynaptic cells of interest to
%            display.  If not given, will default to view all presynaptic cells
%            detected.
% badcnxns (optional) - a list of bad connections listed as pre-post pairs,
%                        one row for each pair, elment 1 of row is pre cell,
%                        element 2 is post.
%
% OUTPUTS
% [userspecified] - outputs a modified version of "funcsynapses" using the
%            name specified in outputvarname
% invisfig - handle to invisible figure that controls the data and function
%            of this entire multi-figure environment
% Brendon Watson 2014

%set up for output
if ~exist('outputvarname','var')
    outputvarname = 'funcsynapses_reviewed';
end
if ~exist('allpresyn','var')%if user did not specify a subset of cells,
    allpresyn = unique(cat(1,funcsynapses.ECells,funcsynapses.ICells));%use all of them
end
if ~exist('badcnxns','var')
    badcnxns = [];
end

% if ~isfield(funcsynapses,'fullRawCCGMtx')
%     funcsynapses = FindSynapse_GetRawCCGMtx(funcsynapses,rawS);
% end

%gather connections
allcnxns = cat(1,funcsynapses.ConnectionsE,funcsynapses.ConnectionsI);

tR = funcsynapses.CCGbins;
allreviewfigs = [];
for a = 1:length(allpresyn);%for each presynaptic cell, we'll make a fig with an axes for each postsyncell
    pre = allpresyn(a);
    thesecnxns = allcnxns(:,1)==pre;
    thesecnxns = allcnxns(thesecnxns,:);
    
%% make figure
    h = figure('Position', get(0,'Screensize'));
    figname = (['ConnectionsPresynCell' num2str(pre)]);
%     set(h,'position',[800 40 1100 925],'name',figname,'CloseRequestFcn',@CloseFigureFcn);
    set(h,'name',figname,'CloseRequestFcn',@CloseFigureFcn);
    allreviewfigs(end+1) = h;
    if ismember(pre,funcsynapses.EIProblemCells)%make pink background as alert if cell has both E and I effects
       set(h,'color',[1 .5 .4])
    end
    
    %get number of subplots per dimension
    subplotdims = ceil(sqrt(size(thesecnxns,1)))+1;

    %% in subplot 1, put autocorrelogram for this cell
    acg = funcsynapses.fullRawCCGMtx(:,pre,pre);
%     T = TimePoints(rawS{1});
%     G = ones(size(T));
%     SampleRate = 10000;
%     BinSize = funcsynapses.BinMs;
%     numBinsBinSize = BinSize*SampleRate/1000;
%     HalfBins = round(300/numBinsBinSize);
%     [acg, dummy, dummy] = CCG(T, G, numBinsBinSize, HalfBins, SampleRate, unique(G), 'count');
    
    sh = subplot(subplotdims,subplotdims,1);
    bar(tR,acg,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])
    axis tight
    
    
    firstplotstring = ['All cnxns where ' num2str(pre) ' is presynaptic']; 
    acgstring = ['Autocorellogram below; ' num2str(funcsynapses.CellRates(pre)) 'Hz'];
    
    if ismember(pre,funcsynapses.EIProblemCells)%if it's a problem cell
            title({firstplotstring,...
                ['THIS PRESYNAPTIC CELL HAS BOTH E AND I CONNECTIONS'],...
                acgstring});
    else%if first subplot but not a problem cell
            title({firstplotstring,...
                acgstring});
    end

    
%% go through each cell with a lagged correlation relative to the "presynaptic" neuron    
    for b = 1:size(thesecnxns,1)
        post = thesecnxns(b,2);
        
        if funcsynapses.FlippedCnxns(pre,post)
            thishibound = funcsynapses.PairUpperThreshs(pre,post);
            thislobound = funcsynapses.PairLowerThreshs(pre,post);
        else
            thishibound = funcsynapses.PairUpperThreshs(post,pre);
            thislobound = funcsynapses.PairLowerThreshs(post,pre);
        end
        
        thisccg = funcsynapses.fullCCGMtx(:,pre,post);
        
% plot ccg
        sh = subplot(subplotdims,subplotdims,b+1);
        bar(tR,thisccg)
        hold on
        
        if funcsynapses.CellShanks(pre) == funcsynapses.CellShanks(post)% if pre & post on same shank...
            % plot correction for spike shadow, which is 0.85ms with process_multi_start
            % 0.85ms is closest 2 spikes can be.  So in binning system with 0.5ms bins, ie from 0.25 to 0.75, that second bin has data    
            % ... just for visualization
            binstarts = tR-funcsynapses.BinMs/2;

            posbin = find(binstarts<0.425,1,'last');
            poscorrection = funcsynapses.BinMs/(0.425-binstarts(posbin));
            posnewval = thisccg(posbin)*poscorrection;
    %         bar(tR(posbin),posnewval,'c')
            plot([tR(posbin) tR(posbin)],[0 posnewval],'LineWidth',1,'color',[0 1 1]);
            plot(tR(posbin), posnewval,'color','c','marker','+');

            negbin = find(binstarts<-0.425,1,'last');
            negcorrection = funcsynapses.BinMs/(binstarts(negbin+1)--0.425);
            negnewval = thisccg(negbin)*negcorrection;
    %         bar(tR(posbin),posnewval,'c')
            plot([tR(negbin) tR(negbin)],[0 negnewval],'LineWidth',1,'color',[0 1 1]);
            plot(tR(negbin), negnewval,'color','c','marker','+');
        end
        
% plot upper and lower significance bounds
        line([min(tR) max(tR)],[thislobound thislobound],'color','m')
        line([min(tR) max(tR)],[thishibound thishibound],'color','m')        
        

% plot E or I star at the proper bin
        if ismember([pre post],funcsynapses.ConnectionsE,'rows')
            thiscolor = 'g';
        elseif ismember([pre post],funcsynapses.ConnectionsI,'rows')
            thiscolor = 'r';
        end
        thispktr = funcsynapses.CnxnTimesVsRefSpk(pre,post);
        thisbin = funcsynapses.CnxnBinsVsCCGStart(pre,post);
        plot(thispktr,thisccg(thisbin),'color',thiscolor,'marker','*','markersize',15)

% plot green and red start and stop vertical lines around the significant bins
        thisstart = funcsynapses.CnxnStartTimesVsRefSpk(pre,post);
        thisstop = funcsynapses.CnxnEndTimesVsRefSpk(pre,post);
        if funcsynapses.FlippedCnxns(pre,post)
            thisstop = -thisstop;
            thisstart = -thisstart;
        end
        binwidth = tR(2)-tR(1);
        plot([thisstart-binwidth/2 thisstart-binwidth/2],[0 thisccg(tR==thisstart)],'color','g','LineWidth',2)
        plot([thisstop+binwidth/2 thisstop+binwidth/2],[0 thisccg(tR==thisstop)],'color','r','LineWidth',1)
        
% plot a vertical line at zero
        xb = get(sh,'XLim');
        yb = get(sh,'YLim');
        plot([0 0],yb,'--c')

% % plot scaled estimate of bins cut by synaptic sampling width
%         if funcsynapses.CellShanks(pre) == funcsynapses.CellShanks(post)
%             1;
%         end
        
%% plot width of any zero-lag correlation detected
%         isE = ismember([pre post],funcsynapses.ZeroLag.EPairs,'rows') || ...
%                 ismember([post pre],funcsynapses.ZeroLag.EPairs,'rows');
%         isI = ismember([pre post],funcsynapses.ZeroLag.IPairs,'rows') || ...
%                 ismember([post pre],funcsynapses.ZeroLag.IPairs,'rows');
%         if isE || isI%bunch of effort to just get the range
%             if isE
%                 [isp,locp] = ismember([pre post],funcsynapses.ZeroLag.EPairs,'rows');
%                 [isn,locn] = ismember([post pre],funcsynapses.ZeroLag.EPairs,'rows');
%                 if isp
%                     thisrange = funcsynapses.ZeroLag.ERanges(locp,:);
%                 elseif isn
%                     thisrange = funcsynapses.ZeroLag.ERanges(locn,:);
%                 end
%             elseif isI
%                 [isp,locp] = ismember([pre post],funcsynapses.ZeroLag.IPairs,'rows');
%                 [isn,locn] = ismember([post pre],funcsynapses.ZeroLag.IPairs,'rows');
%                 if isp
%                     thisrange = funcsynapses.ZeroLag.IRanges(locp,:);
%                 elseif isn
%                     thisrange = funcsynapses.ZeroLag.IRanges(locn,:);
%                 end
%             end
%             if thisrange ~= [0 0];
%                 convertor = funcsynapses.ZeroLag.BinMs/funcsynapses.BinMs;%convert time bins used for zerolag to those for the rest of analysis
%                 thisrange = thisrange*convertor;
% 
%     %             xval = tR(round(mean(thisrange)));
%                 xval = tR(round(thisrange));
%                 if funcsynapses.FlippedCnxns(pre,post)
%                    xval = -fliplr(xval); 
%                 end
%     %             yval = mean(thisccg(round(thisrange)));
%     %             yval = thislobound;
%                 yval = mean(thisccg);
% 
%                 line(xval,[yval yval],'color',[.5 0 1],'LineWidth',5)
%     %             hh = herrorbar(xval,yval,diff(thisrange)/2);
%     %             set(hh,'color',[.5 0 1],'LineWidth',5)
%             end
%         end
%% put title for each cell including cell num, shank and clu num for pre and
% post
        cnxnstring = {['Cell#' num2str(pre) ...
            ' (shank' num2str(funcsynapses.CellShanks(pre)) ...
            ':clu' num2str(funcsynapses.CellShankIDs(pre)) ...
            ') >>']...
            ['Cell#' num2str(post) ...
            ' (shank' num2str(funcsynapses.CellShanks(post)) ...
            ':clu' num2str(funcsynapses.CellShankIDs(post)) ')']};
        title(cnxnstring)

% set up for later management of bad cells or wide-correlation cells        
        bth = text((xb(1)+.15*diff(xb)),yb(1)+.9*(diff(yb)),'','Parent',sh);%blank text to be changed later
        wth = text((xb(1)+.15*diff(xb)),yb(1)+.8*(diff(yb)),'','Parent',sh);%blank text to be changed later
        rh = rectangle('Position',[xb(1) yb(1) diff(xb) diff(yb)],'EdgeColor',[.4 .6 1],'LineWidth',10,'Visible','Off');
        
% store data in userdata for later acces with clicks
        axdata.axes = sh;
        axdata.prepost = [pre post];
        axdata.btexthandle = bth;
        axdata.wtexthandle = wth;
        axdata.recthandle = rh;
        set(sh,'XLim',[min(tR),max(tR)],'UserData',axdata,'ButtonDownFcn',@subplotclick)
        
% Put in an inset figure to display ACG of the postsyn cell
        thisacg = funcsynapses.fullRawCCGMtx(:,post,post);
        axh = AxesInsetBars(sh,.25,[.5 .5 .5],tR',thisacg);
        axhpos = get(axh,'Position');
        set(axh,'Position',[axhpos(1) axhpos(2)-axhpos(4)*.2 axhpos(3) axhpos(4)])
        xb = get(axh,'XLim');
        yb = get(axh,'YLim');
        thisrate = num2str(funcsynapses.CellRates(post));
        if length(thisrate)>4
            thisrate = thisrate(1:4);
        end
        text(xb(1)+0.2*(xb(2)-xb(1)),yb(1)+0.88*(yb(2)-yb(1)),[thisrate 'Hz'])

        if ~isempty(badcnxns)
            if ismember([pre post],badcnxns,'rows')
                set(axdata.axes,'Color',[1 .75 .75])
                set(axdata.btexthandle,'String','BAD')
            end
        end
    end
end

%% Invisible figure for gui management
invisfig = figure('Visible','Off');
for a = allreviewfigs
    set(a,'UserData',invisfig);
end

invisfigdata.allreviewfigs = allreviewfigs;
invisfigdata.badcnxns = badcnxns;
invisfigdata.widecnxns = [];
invisfigdata.outputvarname = outputvarname;
invisfigdata.funcsynapses = funcsynapses;
invisfigdata.thispath = cd;
set(invisfig,'UserData',invisfigdata);




function subplotclick(obj,ev) %when an axes is clicked

figobj = get(obj,'Parent');
invisfig = get(figobj,'UserData');

axdata = get(obj,'UserData');    
prepostcells = axdata.prepost;
btexthandle = axdata.btexthandle;
wtexthandle = axdata.wtexthandle;

invisfigdata = get(invisfig,'UserData');
badcnxns = invisfigdata.badcnxns;

if isfield (invisfigdata,'widecnxns')
    widecnxns = invisfigdata.widecnxns;
else
    widecnxns = [];
end

if strcmp(get(figobj,'SelectionType'),'normal')%if left click, toggle pink/bad pair vs white/ok pair
    clr = get(obj,'Color');
    if sum(clr == [1 1 1])==3%if white (ie synapse), set to pink (bad), remember as bad
        set(obj,'Color',[1 .75 .75])
        badcnxns(end+1,:) = prepostcells;
        set(btexthandle,'String','BAD')
    elseif sum(clr == [1 .75 .75])==3%if pink, set to blue, remember as wide-correlation, unremember as bad
        set(obj,'Color',[1 1 1])
        rowtodelete = find(ismember(badcnxns,prepostcells,'rows'));
    %     widecnxns(rowtodelete,:) = [];    
        badcnxns(rowtodelete,:) = [];    
        set(btexthandle,'String','')
    end
elseif strcmp(get(figobj,'SelectionType'),'alt')%if right click, toggle blue border/wide assn pair vs white/no wide assn pair
    rh = axdata.recthandle;
    switch get(rh,'Visible')
        case 'off'
            set(rh,'Visible','on')
            widecnxns(end+1,:) = prepostcells;
            set(wtexthandle,'String','WIDE')
        case 'on'
            set(rh,'Visible','off')
            rowtodelete = find(ismember(widecnxns,prepostcells,'rows'));
            widecnxns(rowtodelete,:) = [];    
            set(wtexthandle,'String','')
    end
end
        %     rectangle('Position',[xb(1) yb(1) diff(xb) diff(yb)],'EdgeColor',[.6 .6 1],'LineWidth',10)

    
invisfigdata.badcnxns = badcnxns;
invisfigdata.widecnxns = widecnxns;
set(invisfig,'UserData',invisfigdata)


function CloseFigureFcn(obj,ev)%when a figure is closed
invisfig = get(obj,'UserData');
invisfigdata = get(invisfig,'UserData');

allotherfigs = invisfigdata.allreviewfigs;
allotherfigs (find(allotherfigs==obj)) = [];

if ~isempty(allotherfigs)
    name = get(obj,'Name');
    set(obj,'CloseRequestFcn','closereq')

    %save fig in folder ConnectionFigs
    thispath = invisfigdata.thispath;
    if ~exist(fullfile(thispath,'ConnectionFigs'),'dir')
        mkdir(fullfile(thispath,'ConnectionFigs'))
    end
    saveas(obj,fullfile(thispath,'ConnectionFigs',name),'fig');

    delete(obj)
    invisfigdata.allreviewfigs = allotherfigs;
    set(invisfig,'UserData',invisfigdata)
else 
    invisfigdata = get(invisfig,'UserData');
    name = get(obj,'Name');
    
    %save fig in folder ConnectionFigs
    thispath = invisfigdata.thispath;
    if ~exist(fullfile(thispath,'ConnectionFigs'),'dir')
        mkdir(fullfile(thispath,'ConnectionFigs'))
    end
    saveas(obj,fullfile(thispath,'ConnectionFigs',name),'fig');

    delete(obj)
    delete(invisfig)
    
    
    
    %% variable function per user here
    funcsynapses = FindSynapseToStruct(invisfigdata.funcsynapses,invisfigdata.badcnxns,invisfigdata.widecnxns);    
    %%
    assignin('caller',invisfigdata.outputvarname,funcsynapses)

    oi = findobj('Type','Figure','Visible','Off','Name','BWWaitFig');
    set(oi,'Name','DELETEMENOW')
end
[];

