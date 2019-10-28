function [] = LogScale( whichaxis,logbase,varargin)
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
p = inputParser;
addParameter(p,'exp',false)
addParameter(p,'nohalf',false)
parse(p,varargin{:})
expo = p.Results.exp;
nohalf = p.Results.nohalf;

%%
if strcmp(whichaxis,'y') || strcmp(whichaxis,'xy')
    range = get(gca,'YLim');
    switch nohalf
        case true
            range(1) = ceil(range(1)); range(2) = floor(range(2));
        case false
            range(1) = ceil(range(1).*2)./2; range(2) = floor(range(2).*2)./2;
    end
    ticks = [range(1):range(2)];
    if length(ticks)<=3 & ~nohalf
        ticks = [range(1):0.5:range(2)];
    end
    
    if length(ticks)>5
        ticks = ticks(1:2:end);
    end
    
    set(gca,'YTick',ticks)
    set(gca,'YTickLabels',round(logbase.^ticks,2,'significant'))
    
    if expo && logbase==10
        tickstrings = cellfun(@num2str,(num2cell(ticks)),'uniformoutput',false);
        tickstrings = cellfun(@(X) replace(X,'','^'),tickstrings,'uniformoutput',false);
        tickstrings = cellfun(@(X) [num2str(logbase),X(1:end-1)],tickstrings,'uniformoutput',false);
        set(gca,'YTickLabels',tickstrings)
    end
    
end

if strcmp(whichaxis,'x') || strcmp(whichaxis,'xy')
    range = get(gca,'XLim');
    switch nohalf
        case true
            range(1) = ceil(range(1)); range(2) = floor(range(2));
        case false
            range(1) = ceil(range(1).*2)./2; range(2) = floor(range(2).*2)./2;
    end
    ticks = [range(1):range(2)];
    if length(ticks)<=3 & ~nohalf
        ticks = [range(1):0.5:range(2)];
    end
    
    if length(ticks)>5
        ticks = ticks(1:2:end);
    end
        
    
    set(gca,'XTick',ticks)
    set(gca,'XTickLabels',round(logbase.^ticks,2,'significant'))
    
    if expo && logbase==10
        tickstrings = cellfun(@num2str,(num2cell(ticks)),'uniformoutput',false);
        tickstrings = cellfun(@(X) replace(X,'','^'),tickstrings,'uniformoutput',false);
        tickstrings = cellfun(@(X) [num2str(logbase),X(1:end-1)],tickstrings,'uniformoutput',false);
        set(gca,'XTickLabels',tickstrings)
    end
end

if strcmp(whichaxis,'z')
    range = get(gca,'ZLim');
    switch nohalf
        case true
            range(1) = ceil(range(1)); range(2) = floor(range(2));
        case false
            range(1) = ceil(range(1).*2)./2; range(2) = floor(range(2).*2)./2;
    end
    ticks = [range(1):range(2)];
    if length(ticks)<=3 & ~nohalf
        ticks = [range(1):0.5:range(2)];
    end
    
    if length(ticks)>5
        ticks = ticks(1:2:end);
    end
    
    set(gca,'ZTick',ticks)
    set(gca,'ZTickLabels',round(logbase.^ticks,2,'significant'))
    
    if expo
        tickstrings = cellfun(@num2str,(num2cell(ticks)),'uniformoutput',false);
        tickstrings = cellfun(@(X) replace(X,'','^'),tickstrings,'uniformoutput',false);
        tickstrings = cellfun(@(X) [num2str(logbase),X(1:end-1)],tickstrings,'uniformoutput',false);
        set(gca,'ZTickLabels',tickstrings)
    end
end


if strcmp(whichaxis,'c')
    range = get(gca,'CLim');
    switch nohalf
        case true
            range(1) = ceil(range(1)); range(2) = floor(range(2));
        case false
            range(1) = ceil(range(1).*2)./2; range(2) = floor(range(2).*2)./2;
    end
    ticks = [range(1):range(2)];
    if length(ticks)<=3 & ~nohalf
        ticks = [range(1):0.5:range(2)];
    end
    
    if length(ticks)>=5
        ticks = ticks(1:2:end);
    end
    
    cb = get(gca,'colorbar');
    
    set(cb,'Ticks',ticks)
    set(cb,'TickLabels',round(logbase.^ticks,1,'significant'))
    
    if expo
        tickstrings = cellfun(@num2str,(num2cell(ticks)),'uniformoutput',false);
        tickstrings = cellfun(@(X) replace(X,'','^'),tickstrings,'uniformoutput',false);
        tickstrings = cellfun(@(X) [num2str(logbase),X(1:end-1)],tickstrings,'uniformoutput',false);
        set(cb,'TickLabels',tickstrings)
    end
    
end

end

