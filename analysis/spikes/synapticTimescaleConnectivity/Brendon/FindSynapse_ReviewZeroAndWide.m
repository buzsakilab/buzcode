function FindSynapse_ReviewZeroAndWide(funcsynapses,outputvarname)
% Takes output of FindSynapse.m and plots positive connections as ccr
% plots for review and clickable deletions of bad connections.
% Brendon Watson 2014

if ~exist('outputvarname','var')
    outputvarname = 'funcsynapses_reviewed';
end

ze = unique(funcsynapses.ZeroLag.EPairs,'rows');
zi = unique(funcsynapses.ZeroLag.IPairs,'rows');

w = unique(funcsynapses.WideConnections,'rows');

% allpresyn = unique(cat(1,funcsynapses.ECells,funcsynapses.ICells));
% allcnxns = cat(1,funcsynapses.ConnectionsE,funcsynapses.ConnectionsI);

tR = funcsynapses.CCGbins;
allreviewfigs = [];
numcells = size(funcsynapses.CellRates,1);
for a = 1:numcells;
    pre = a;
    
    if ismember(pre,ze) | ismember(pre,zi) | ismember(pre,w)
        [ix1,ix2] = find(ze == pre);
        theseze = ze(ix1,:);
        
        [ix1,ix2] = find(zi == pre);
        thesezi = zi(ix1,:);
        
        [ix1,ix2] = find(w == pre);
        thesew = w(ix1,:);

        thesecnxns = cat(1,theseze,thesezi,thesew);
        eiw = cat(1,1*ones(size(theseze,1),1),2*ones(size(thesezi,1),1),3*ones(size(thesew,1),1));
        
        h = figure;
        figname = (['ConnectionsPresynCell' num2str(pre)]);
        set(h,'position',[800 40 1100 925],'name',figname);
%         set(h,'CloseRequestFcn',@CloseFigureFcn);
        allreviewfigs(end+1) = h;
        if ismember(pre,funcsynapses.EIProblemCells)
           set(h,'color',[1 .5 .4])
        end

        subplotdims = ceil(sqrt(size(thesecnxns,1)))+1;

        acg = funcsynapses.fullRawCCGMtx(:,pre,pre);
        sh = subplot(subplotdims,subplotdims,1);
        bar(tR,acg,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])
        axis tight


        firstplotstring = ['All Zero or Wide where ' num2str(pre) ' is involved']; 
        acgstring = ['Autocorellogram below; ' num2str(funcsynapses.CellRates(pre)) 'Hz'];

%         if ismember(pre,funcsynapses.EIProblemCells)%if it's a problem cell
%                 title({firstplotstring,...
%                     ['THIS PRESYNAPTIC CELL HAS BOTH E AND I CONNECTIONS'],...
%                     acgstring});
%         else%if first subplot but not a problem cell
                title({firstplotstring,...
                    acgstring});
%         end



        for b = 1:size(thesecnxns,1)
%             post = thesecnxns(b,2);
            ix = thesecnxns(b,:)~=pre;
            post = thesecnxns(b,ix);

            if funcsynapses.FlippedCnxns(pre,post)
                thishibound = funcsynapses.PairUpperThreshs(pre,post);
                thislobound = funcsynapses.PairLowerThreshs(pre,post);
            else
                thishibound = funcsynapses.PairUpperThreshs(post,pre);
                thislobound = funcsynapses.PairLowerThreshs(post,pre);
            end

            thisccg = funcsynapses.fullCCGMtx(:,pre,post);

            sh = subplot(subplotdims,subplotdims,b+1);
            bar(tR,thisccg)
            hold on
            
            if funcsynapses.CellShanks(pre) == funcsynapses.CellShanks(post)% if pre & post on same shank
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
            
            line([min(tR) max(tR)],[thislobound thislobound],'color','m')
            line([min(tR) max(tR)],[thishibound thishibound],'color','m')        

            thispktr = funcsynapses.CnxnTimesVsRefSpk(pre,post);
            thisbin = funcsynapses.CnxnBinsVsCCGStart(pre,post);

            if ~isempty(funcsynapses.ConnectionsE)
                if ismember([pre post],funcsynapses.ConnectionsE,'rows')
                    thiscolor = 'g';
                end
            elseif ~isempty(funcsynapses.ConnectionsI)
                if ismember([pre post],funcsynapses.ConnectionsI,'rows')
                    thiscolor = 'r';
                end
            end
%             plot(thispktr,thisccg(thisbin),'color',thiscolor,'marker','*','markersize',15)

            xb = get(sh,'XLim');
            yb = get(sh,'YLim');
            plot([0 0],yb,'--c')

            cnxnstring = {['Cell#' num2str(pre) ...
                ' (clu' num2str(funcsynapses.CellShankIDs(pre)) ...
                ':shank' num2str(funcsynapses.CellShanks(pre)) ...
                ') >>']...
                ['Cell#' num2str(post) ...
                ' (clu' num2str(funcsynapses.CellShankIDs(post)) ...
                ':shank' num2str(funcsynapses.CellShanks(post)) ')']};
            title(cnxnstring)
            th = text((xb(1)+.15*diff(xb)),yb(1)+.9*(diff(yb)),'','Parent',sh);%blank text to be changed later

%             axdata.prepost = [pre post];
%             axdata.texthandle = th;
%             set(sh,'UserData',axdata,'ButtonDownFcn',@subplotclick)
            set(sh,'XLim',[min(tR),max(tR)])

            %% Put in an inset figure
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
            
            %% Make clear whether was zerolag or wide
            switch eiw(b)
                case(1)
                    set(th,'String','ZeroE')
                case(2)
                    set(th,'String','ZeroI')
                case(3)
                    set(th,'String','Wide')
            end            
        end
        
        thispath = cd;
        if ~exist(fullfile(thispath,'ZeroLagAndWideFigs'),'dir')
            mkdir(fullfile(thispath,'ZeroLagAndWideFigs'))
        end
        saveas(h,fullfile(thispath,'ZeroLagAndWideFigs',figname),'fig');
    end
end

% h1 = figure('Visible','Off');
% for a = allreviewfigs
%     set(a,'UserData',h1);
% end

% h1data.allreviewfigs = allreviewfigs;
% h1data.badcnxns = [];
% h1data.widecnxns = [];
% h1data.outputvarname = outputvarname;
% h1data.funcsynapses = funcsynapses;
% set(h1,'UserData',h1data);


% 
% 
% function subplotclick(obj,ev)
% 
% figobj = get(obj,'Parent');
% invisfig = get(figobj,'UserData');
% 
% axdata = get(obj,'UserData');    
% prepostcells = axdata.prepost;
% texthandle = axdata.texthandle;
% 
% h1data = get(invisfig,'UserData');
% badcnxns = h1data.badcnxns;
% 
% if isfield (h1data,'widecnxns')
%     widecnxns = h1data.widecnxns;
% else
%     widecnxns = [];
% end
% 
% clr = get(obj,'Color');
% if sum(clr == [1 1 1])==3%if white (ie synapse), set to pink (bad), remember as bad
%     set(obj,'Color',[1 .75 .75])
%     badcnxns(end+1,:) = prepostcells;
%     set(texthandle,'String','BAD')
% elseif sum(clr == [1 .75 .75])==3%if pink, set to blue, remember as wide-correlation, unremember as bad
%     set(obj,'Color',[.75 .75 1])
%     widecnxns(end+1,:) = prepostcells;
%     rowtodelete = find(ismember(badcnxns,prepostcells,'rows'));
%     badcnxns(rowtodelete,:) = [];     
%     set(texthandle,'String','WIDE')
% elseif sum(clr == [.75 .75 1])==3%if is blue, set to wite, unremember as wide-correlation
%     set(obj,'Color',[1 1 1])
%     rowtodelete = find(ismember(badcnxns,prepostcells,'rows'));
%     widecnxns(rowtodelete,:) = [];    
%     set(texthandle,'String','')
% end
% 
% h1data.badcnxns = badcnxns;
% h1data.widecnxns = widecnxns;
% set(invisfig,'UserData',h1data)
% 
% 
% function CloseFigureFcn(obj,ev)
% invisfig = get(obj,'UserData');
% h1data = get(invisfig,'UserData');
% 
% allotherfigs = h1data.allreviewfigs;
% allotherfigs (find(allotherfigs==obj)) = [];
% 
% if ~isempty(allotherfigs)
%     delete(obj)
%     h1data.allreviewfigs = allotherfigs;
%     set(invisfig,'UserData',h1data)
% else 
%     h1data = get(invisfig,'UserData');
%     delete(obj)
%     delete(invisfig)
%     
%     %% variable function per user here
%     funcsynapses = FindSynapseToStruct(h1data.funcsynapses,h1data.badcnxns,h1data.widecnxns);    
%     %%
%     assignin('caller',h1data.outputvarname,funcsynapses)
% end
% [];

