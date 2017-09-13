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
    range(1) = round(range(1)); range(2) = round(range(2));
    ticks = [range(1):range(2)];
    if length(ticks)<=4
        ticks = [range(1):0.5:range(2)];
    end
    
    set(gca,'YTick',ticks)
    set(gca,'YTickLabels',round(logbase.^ticks,3,'significant'))
end

if strcmp(whichaxis,'x') || strcmp(whichaxis,'xy')
    range = get(gca,'XLim');
    range(1) = round(range(1)); range(2) = round(range(2));
    ticks = [range(1):range(2)];
    if length(ticks)<=4
        ticks = [range(1):0.5:range(2)];
    end
    set(gca,'XTick',ticks)
    set(gca,'XTickLabels',round(logbase.^ticks,3,'significant'))
end

if strcmp(whichaxis,'z')
    range = get(gca,'ZLim');
    range(1) = round(range(1)); range(2) = round(range(2));
    ticks = [range(1):range(2)];
    if length(ticks)<=4
        ticks = [range(1):0.5:range(2)];
    end
    set(gca,'ZTick',ticks)
    set(gca,'ZTickLabels',round(logbase.^ticks,3,'significant'))
end

end

