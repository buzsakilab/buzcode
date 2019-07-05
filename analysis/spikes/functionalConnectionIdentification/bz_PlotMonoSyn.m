function keep_con = bz_PlotMonoSyn(ccgR,sig_con,Pred,Bounds,completeIndex,binSize,duration)
% Manual sorting
% click to deselect the ccg (turns pink)

% Edited by Peter Petersen

keep_con = sig_con;
window  =false(size(ccgR,1),1);
window(ceil(length(window)/2) - round(.004/binSize): ceil(length(window)/2) + round(.004/binSize)) = true;
halfBins = round(duration/binSize/2);

t = (-halfBins:halfBins)'*binSize;

allcel = unique(sig_con(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure('KeyReleaseFcn', {@keyPress},'Name','MonoSynCon inspector','NumberTitle','off','renderer','opengl');

for i=1:length(allcel)
    if ~ishandle(h)
        return
    end
    cla, clf
    
    prs = sig_con(any(sig_con==allcel(i),2),:);
    ha = tight_subplot(ceil(sqrt(1+size(prs,1))),ceil(sqrt(1+size(prs,1))),[.02 .01]);
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
            bar_from_patch(t,ccgR(:,prs1(1),prs1(2)),'b')
            % bar(t,ccgR(:,prs1(1),prs1(2)),1,'FaceColor','b','EdgeColor','b');
            hold on;
            
            % Plot predicted values
            plot(t,Pred(:,prs1(1),prs1(2)),'g');
            
            %Plot upper and lower boundaries
            plot(t,Bounds(:,prs1(1),prs1(2),1),'r--');
            plot(t,Bounds(:,prs1(1),prs1(2),2),'r--');
            
            % Plot signif increased bins in red
            bar_from_patch(t,exc,'r')
            % bar(t,exc,1,'FaceColor','r','EdgeColor','r');
            
            % Plot signif lower bins in blue
            bar_from_patch(t,inh,'c')
            % bar(t,inh,1,'FaceColor','c','EdgeColor','c');
            
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
            bar_from_patch(t,ccgR(:,allcel(i),allcel(i)),'k')
            % bar(t,ccgR(:,allcel(i),allcel(i)),1,'FaceColor','k','EdgeColor','k');
            xlim([min(t) max(t)]);
            xlabel('Reference Cell ACG');
            targ=completeIndex(completeIndex(:,3) == allcel(i),1:2);
            mtit(['Reference Cell sh: ' num2str(targ(1)) ' cell '  num2str(targ(2)),' (', num2str(i),'/' num2str(length(allcel)),')'])
            
            uiwait(h);
        end
    end
    
end
if ishandle(h)
    close(h)
end

    function subplotclick(obj,ev) %when an axes is clicked
        figobj = get(obj,'Parent');
        axdata = get(obj,'UserData');
        clr2 = get(obj,'Color');
        if sum(clr2 == [1 1 1])==3%if white (ie synapse), set to pink (bad), remember as bad
            set(obj,'Color',[1 .75 .75])
            keep_con(ismember(keep_con,prs(axdata,:),'rows'),:)=[];
        elseif sum(clr2 == [1 .75 .75])==3%if pink, set to white, set to good
            set(obj,'Color',[1 1 1])
            keep_con = [keep_con;prs(axdata,:)];
        end
    end

    function advance
        uiresume(h);
    end

    function keyPress(src, e)
        switch e.Key
            case 'space'
                uiresume(h);
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

bar_from_patch(xdata,ydata,color)
% bar(xdata,ydata,1,'FaceColor',color,'EdgeColor',color)
axis tight
end

function bar_from_patch(x_data, y_data,col)
% Creates a bar graph using the patch plot mode, which is substantial
% faster than using the regular bar plot. 
% By Peter Petersen

x_step = x_data(2)-x_data(1);
x_data = [x_data(1),reshape([x_data,x_data+x_step]',1,[]),x_data(end)];
y_data = [0,reshape([y_data,y_data]',1,[]),0];
patch(x_data, y_data,col,'EdgeColor',col)
end