function [] = LogScale_ss( whichaxis,logbase)
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
    range(1) = floor(range(1));range(2) = ceil(range(2));
    set(gca,'YTick',[range(1):range(2)])
    set(gca,'YTickLabels',logbase.^[range(1):range(2)])
end

if strcmp(whichaxis,'x') || strcmp(whichaxis,'xy')
    range = get(gca,'XLim');
    range(1) = floor(range(1));range(2) = ceil(range(2));
    set(gca,'XTick',[range(1):range(2)])
    set(gca,'XTickLabels',logbase.^[range(1):range(2)])
end

if strcmp(whichaxis,'z')
    range = get(gca,'ZLim');
    range(1) = floor(range(1));range(2) = ceil(range(2));
    set(gca,'ZTick',[range(1):range(2)])
    set(gca,'ZTickLabels',logbase.^[range(1):range(2)])
end

end

