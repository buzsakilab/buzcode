function [] = LogScale( whichaxis,logbase)
%LogScale(whichaxis,logbase) renames the tick labels on your axes for a
%logarithmically scaled variable.
%
%INPUT
%   whichaxis   which axis did you log scale (can be 'x' 'y' or 'xy')
%   logbase     base of the logarithm you used for your points (ex 10)
%
%EXAMPLE
%   plot(log10(Xvals),Yvals)
%   LogScale('x',10)
%
%DLevenstein 2015
%%

if strcmp(whichaxis,'y') || strcmp(whichaxis,'xy')
    range = get(gca,'YLim');
    range(1) = ceil(range(1)); range(2) = floor(range(2));
    ticks = [range(1):range(2)];
    if length(ticks)<=3
        ticks = [range(1):0.5:range(2)];
    end
    
    if length(ticks)>=5
        ticks = ticks(1:2:end);
    end
    
    set(gca,'YTick',ticks)
    set(gca,'YTickLabels',round(logbase.^ticks,3,'significant'))
    
    if max(abs(ticks))>=3
        tickstrings = cellfun(@num2str,(num2cell(ticks)),'uniformoutput',false);
        tickstrings = cellfun(@(X) replace(X,'','^'),tickstrings,'uniformoutput',false);
        tickstrings = cellfun(@(X) ['10',X(1:end-1)],tickstrings,'uniformoutput',false);
        set(gca,'YTickLabels',tickstrings)
    end
    
end

if strcmp(whichaxis,'x') || strcmp(whichaxis,'xy')
    range = get(gca,'XLim');
    range(1) = ceil(range(1)); range(2) = floor(range(2));
    ticks = [range(1):range(2)];
    if length(ticks)<=3
        ticks = [range(1):0.5:range(2)];
    end
    
    if length(ticks)>=5
        ticks = ticks(1:2:end);
    end
        
    
    set(gca,'XTick',ticks)
    set(gca,'XTickLabels',round(logbase.^ticks,3,'significant'))
    
    if max(abs(ticks))>=3
        tickstrings = cellfun(@num2str,(num2cell(ticks)),'uniformoutput',false);
        tickstrings = cellfun(@(X) replace(X,'','^'),tickstrings,'uniformoutput',false);
        tickstrings = cellfun(@(X) ['10',X(1:end-1)],tickstrings,'uniformoutput',false);
        set(gca,'XTickLabels',tickstrings)
    end
end

if strcmp(whichaxis,'z')
    range = get(gca,'ZLim');
    range(1) = ceil(range(1)); range(2) = floor(range(2));
    ticks = [range(1):range(2)];
    if length(ticks)<=3
        ticks = [range(1):0.5:range(2)];
    end
    
    if length(ticks)>=5
        ticks = ticks(1:2:end);
    end
    
    set(gca,'ZTick',ticks)
    set(gca,'ZTickLabels',round(logbase.^ticks,3,'significant'))
    
    if max(abs(ticks))>=2
        tickstrings = cellfun(@num2str,(num2cell(ticks)),'uniformoutput',false);
        tickstrings = cellfun(@(X) replace(X,'','^'),tickstrings,'uniformoutput',false);
        tickstrings = cellfun(@(X) ['10',X(1:end-1)],tickstrings,'uniformoutput',false);
        set(gca,'ZTickLabels',tickstrings)
    end
end


if strcmp(whichaxis,'c')
    range = get(gca,'CLim');
    range(1) = ceil(range(1)); range(2) = floor(range(2));
    ticks = [range(1):range(2)];
    if length(ticks)<=3
        ticks = [range(1):0.5:range(2)];
    end
    
    if length(ticks)>=5
        ticks = ticks(1:2:end);
    end
    
    cb = get(gca,'colorbar');
    
    set(cb,'Ticks',ticks)
    set(cb,'TickLabels',round(logbase.^ticks,3,'significant'))
    
    if max(abs(ticks))>=2
        tickstrings = cellfun(@num2str,(num2cell(ticks)),'uniformoutput',false);
        tickstrings = cellfun(@(X) replace(X,'','^'),tickstrings,'uniformoutput',false);
        tickstrings = cellfun(@(X) ['10',X(1:end-1)],tickstrings,'uniformoutput',false);
        set(cb,'TickLabels',tickstrings)
    end
    
end

end

