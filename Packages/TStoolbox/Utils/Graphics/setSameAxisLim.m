function setSameAxisLim(axs,  directions)
% sets the same axis limits on several plots
%
% INPUTS:
% axs: a cell array of axis handles
% directions: a string composed of 'x' and 'y', indicating which directions
% to act on
% axes will be set to  the maximum interval including all axis current intervals 



    xa = zeros(0,2);
    ya = zeros(0,2);
    for i = 1:length(axs)
        ax = axs{i};
        xa = [xa ; get(ax, 'xlim')];
        ya = [ya ; get(ax, 'ylim')];
    end
    
    minx = min(xa);
    maxx = max(xa);
    intx = [minx(1), maxx(2)];
    
    miny = min(ya);
    maxy = max(ya);
    inty = [miny(1) maxy(2)];
    
    for i = 1:length(axs)
        ax = axs{i};
        if ismember('x', directions)
            set(ax, 'xlim', intx);
        end
        if ismember('y', directions)
            set(ax, 'ylim', inty);
        end
    end
    