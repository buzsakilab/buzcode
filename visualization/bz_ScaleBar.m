function [  ] = bz_ScaleBar( units )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        hold on
       % units = 'AU';
        yrange = get(gca,'ylim');
        xrange = get(gca,'xlim');
        xticks = get(gca,'xtick');
        xbar = xticks(1:2);xbarsize = diff(xbar);
        set(gca,'xtick',xrange-0.5.*xbarsize)
        set(gca,'xticklabels',[num2str(xbarsize),' ',units])
        plot(xrange,[1 1].*yrange(1),'-w','linewidth',3)
        plot(xrange(2)-[xbarsize 0],[1 1].*yrange(1),'-k','linewidth',3)
end

