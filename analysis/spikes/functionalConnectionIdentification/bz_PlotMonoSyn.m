function keep_con = bz_PlotMonoSyn(ccgR,sig_con,Pred,Bounds,completeIndex,binSize,duration)
% Manual curating detected CCGs
% click to deselect the ccg (turns pink)

% Edited by Peter Petersen
% petersen.peter@gmail.com
% Last edited: 12-07-2019

keep_con = sig_con;
window  =false(size(ccgR,1),1);
window(ceil(length(window)/2) - round(.004/binSize): ceil(length(window)/2) + round(.004/binSize)) = true;
halfBins = round(duration/binSize/2);

t = (-halfBins:halfBins)'*binSize;

allcel = unique(sig_con(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure('KeyReleaseFcn', {@keyPress},'Name','MonoSynCon inspector','NumberTitle','off','renderer','opengl');
% p = uipanel(h,'Position',[0 0 1 .1],'BorderType','none')
% p2 = uipanel(h,'Position',[0 0 0.01 0.01],'BorderType','none')
uicontrol('Parent',h,'Style','pushbutton','Position',[5 405 20 15],'Units','normalized','String','<','Callback',@(src,evnt)goBack,'KeyPressFcn', {@keyPress});
uicontrol('Parent',h,'Style','pushbutton','Position',[540 405 20 15],'Units','normalized','String','>','Callback',@(src,evnt)advance,'KeyPressFcn', {@keyPress});
plotTitle = uicontrol('Parent',h,'Style','text','Position',[130 410 350 10],'Units','normalized','String','','HorizontalAlignment','center','FontSize',13);

i = 1;
while i > 0 && i <= length(allcel)
    if ~ishandle(h)
        return
    end
    delete(findobj(h, 'type', 'axes'));
    
    prs = sig_con(any(sig_con==allcel(i),2),:);
    [plotRows,~]= numSubplots(max(2+size(prs,1),4));
    ha = tight_subplot(plotRows(1),plotRows(2),[.02 .01]);
    prs2 = [];
    for j=1:length(ha)
        axes(ha(j))
        if j<=size(prs,1)
            prs1 = prs(j,:);
            if prs1(1)~=allcel(i)
                prs1 = fliplr(prs1);
            end
            prs2(j,:) = prs1;
            exc=ccgR(:,prs1(1),prs1(2));
            exc(exc<Bounds(:,prs1(1),prs1(2),1)|~window)=0;
            
            inh=ccgR(:,prs1(1),prs1(2));
            inh(inh>Bounds(:,prs1(1),prs1(2),2)|~window)=0;
            bar_from_patch(t,ccgR(:,prs1(1),prs1(2)),'b')
            % bar(t,ccgR(:,prs1(1),prs1(2)),1,'FaceColor','b','EdgeColor','b');
            hold on;
            
            % Plot predicted values
            plot(t,Pred(:,prs1(1),prs1(2)),'g', 'HitTest','off');
            
            %Plot upper and lower boundaries
            plot(t,Bounds(:,prs1(1),prs1(2),1),'r--', 'HitTest','off');
            plot(t,Bounds(:,prs1(1),prs1(2),2),'r--', 'HitTest','off');
            
            % Plot signif increased bins in red
            bar_from_patch(t,exc,'r')
            % bar(t,exc,1,'FaceColor','r','EdgeColor','r');
            
            % Plot signif lower bins in blue
            bar_from_patch(t,inh,'c')
            % bar(t,inh,1,'FaceColor','c','EdgeColor','c');
            
            upL = get(gca,'ylim');
            plot([0 0],[0 upL(2)],'k', 'HitTest','off')
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
            
        elseif j<length(ha)-1
            axis off
        elseif j<length(ha)
            zdata = ccgR(:,:,allcel(i))';
            imagesc(flip(t),1:size(ccgR,3),(zdata'./max(zdata'))'), hold on
            plot(-0.058*ones(size(prs2,1),1),prs2(:,2),'.w', 'HitTest','off', 'MarkerSize',12)
            plot(-0.058*ones(size(prs2,1),1),prs2(:,2),'ok', 'HitTest','off')
            plot(-0.058*ones(size(prs2,1),1),prs2(:,1),'.k', 'HitTest','off', 'MarkerSize',12)
            plot(-0.058*ones(size(prs2,1),1),prs2(:,1),'ok', 'HitTest','off')
        else
            bar_from_patch(t,ccgR(:,allcel(i),allcel(i)),'k')
            % bar(t,ccgR(:,allcel(i),allcel(i)),1,'FaceColor','k','EdgeColor','k');
            xlim([min(t) max(t)]);
            xlabel('Reference Cell ACG');
            targ=completeIndex(completeIndex(:,3) == allcel(i),1:2);
            plotTitle.String = ['Reference Cell sh: ' num2str(targ(1)) ' cell '  num2str(targ(2)),' (', num2str(i),'/' num2str(length(allcel)),')'];
            %             mtit(['Reference Cell sh: ' num2str(targ(1)) ' cell '  num2str(targ(2)),' (', num2str(i),'/' num2str(length(allcel)),')'])
            
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

    function goBack
        i = max(i-1,1);
        uiresume(h);
    end

    function advance
        if i==length(allcel)
            answer = questdlg('All cells have been currated. Do you want to quit?', 'Monosyn curration complete', 'Yes','No','Yes');
            if strcmp(answer,'Yes')
                close(h)
            end
        else
            i = i+1;
            uiresume(h);
        end
    end

    function keyPress(src, e)
        switch e.Key
            case 'space'
                advance
            case 'rightarrow'
                advance
            case 'leftarrow'
                goBack
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
axh = axes('Position',newpos, 'HitTest','off');

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
patch(x_data, y_data,col,'EdgeColor',col, 'HitTest','off')
end

function [p,n]=numSubplots(n)
% Calculate how many rows and columns of sub-plots are needed to
% neatly display n subplots.
% Rob Campbell - January 2010

while isprime(n) && n>4
    n=n+1;
end
p=factor(n);
if length(p)==1
    p=[1,p];
    return
end
while length(p)>2
    if length(p)>=4
        p(1)=p(1)*p(end-1);
        p(2)=p(2)*p(end);
        p(end-1:end)=[];
    else
        p(1)=p(1)*p(2);
        p(2)=[];
    end
    p=sort(p);
end

%Reformat if the column/row ratio is too large: we want a roughly
%square design
while p(2)/p(1)>2.5
    N=n+1;
    [p,n]=numSubplots(N); %Recursive!
end
end