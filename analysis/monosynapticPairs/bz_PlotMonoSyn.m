function keep_con = bz_PlotMonoSyn(ccgR,sig_con,Pred,Bounds,completeIndex,binSize,duration)


% Manual sorting
%click to delete the ccg (turns pink)


keep_con = sig_con;
window  =false(size(ccgR,1),1);
window(ceil(length(window)/2) - round(.004/binSize): ceil(length(window)/2) + round(.004/binSize)) = true;
halfBins = round(duration/binSize/2);

t = (-halfBins:halfBins)'*binSize;

allcel = unique(sig_con(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(allcel)
    prs = sig_con(any(sig_con==allcel(i),2),:);
    
    
    h = figure;
    ha = tight_subplot(ceil(sqrt(1+size(prs,1))),ceil(sqrt(1+size(prs,1))),[.02 .01]);
    set(h,'position',[ 1           1        1920        1001])
    
    for j=1:length(ha)
           axes(ha(j))
        if j<=size(prs,1)
            prs1 = prs(j,:);
            if prs1(1)~=allcel(i)
                prs1 = fliplr(prs1);
            end
            
         
            
            
            exc=ccgR(:,prs1(1),prs1(2));
            exc(exc<Bounds(:,prs1(1),prs1(2),1)|~window)=0;
            
            inh=ccgR(:,prs1(1),prs1(2));
            inh(inh>Bounds(:,prs1(1),prs1(2),2)|~window)=0;
            
            
            
            bar(t,ccgR(:,prs1(1),prs1(2)));
            hold on;
            % Plot predicted values
            plot(t,Pred(:,prs1(1),prs1(2)),'g');
            %Plot upper and lower boundaries
            plot(t,Bounds(:,prs1(1),prs1(2),1),'r--');
            plot(t,Bounds(:,prs1(1),prs1(2),2),'r--');
            
            
            % Plot signif increased bins in red
            
            
            
            bar(t,exc,'r');
            % Plot signif lower bins in blue
            
            bar(t,inh,'c');
            
            
            
            
            
            upL = get(gca,'ylim');
            plot([0 0],[0 upL(2)],'k')
            xlim([min(t) max(t)]);
            
            text(min(t) +.1*abs(min(t)),double(max(ccgR(:,prs(j,1),prs(j,2)))),['max cnt: ' num2str(max(ccgR(:,prs(j,1),prs(j,2))))])
            
            
            set(gca,'yticklabel',[],'xtick',[min(t) 0 max(t)],'xticklabel',[])
            
            tcel = setdiff(prs(j,:),allcel(i,:));
            targ=completeIndex(completeIndex(:,3)==tcel,1:2);
            xlabel(['sh: ' num2str(targ(1)) ' cell '  num2str(targ(2))]);
            
            
            %the bad ones are in pink
            if  ~ismember(prs(j,:),keep_con,'rows')
                set(ha(j),'Color',[1 .75 .75])
                
            end
            
            set(gca,'UserData',j,'ButtonDownFcn',@subplotclick);
            
            
            
            % Plot an inset with the ACG
            
            
            thisacg = ccgR(:,tcel,tcel);
            
            
            axh = AxesInsetBars(gca,.2,[.5 .5 .5],t,thisacg);
            axhpos = get(axh,'Position');
            set(axh,'Position',[axhpos(1) axhpos(2)-axhpos(4)*.2 axhpos(3) axhpos(4)],'XTickLabel',[]);
            
            
        elseif j<length(ha)
            axis off
            
        else
            
            
            
            bar(t,ccgR(:,allcel(i),allcel(i)),'k');
            xlim([min(t) max(t)]);
            xlabel('Reference Cell ACG');
            
            
            targ=completeIndex(completeIndex(:,3) == allcel(i),1:2);
            
            
            mtit(['Reference Cell sh: ' num2str(targ(1)) ' cell '  num2str(targ(2))])
            
            
            waitfor(h)
            
        end
    end
    
end


    function subplotclick(obj,ev) %when an axes is clicked
        
        figobj = get(obj,'Parent');
        axdata = get(obj,'UserData');
        
        
        
        clr = get(obj,'Color');
        if sum(clr == [1 1 1])==3%if white (ie synapse), set to pink (bad), remember as bad
            set(obj,'Color',[1 .75 .75])
            
            
            
            
          
            
               keep_con(ismember(keep_con,prs(axdata,:),'rows'),:)=[];
            
            
        elseif sum(clr == [1 .75 .75])==3%if pink, set to white, set to good
            set(obj,'Color',[1 1 1])
           keep_con = [keep_con;prs(axdata,:)];
        end
    end
end


function axh = AxesInsetBars(h,ratio,color,xdata,ydata)
% function axh = AxesInsetBars(h,ratio,color,xdata,ydata)
% Puts an inset bar plot (axes) into the upper right of the given axes.
%  INPUTS
%  h = handle of reference axes
%  ratio = Size of inset relative to original plot
%  color = color of bars
%  data = data to plot
%
%  OUTPUTS
%  axh = handle of inset axes


figpos = get(h,'Position');
newpos = [figpos(1)+(1-ratio)*figpos(3) figpos(2)+(1-ratio)*figpos(4) ratio*figpos(3) ratio*figpos(4)];
axh = axes('Position',newpos);

bar(xdata,ydata,'FaceColor',color,'EdgeColor',color)
axis tight
end
