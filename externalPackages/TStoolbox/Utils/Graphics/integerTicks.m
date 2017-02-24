function integerTicks(ax, directions)
% set x/yticks to  integer values only
%
% INPUTS
% ax: an axis handle
% directions: a string containing 'x' and 'y' indicating which directions
% to act on 
if ismember('x', directions)
    set(ax, 'XTickMode', 'auto');
    yl = get(ax, 'XTick');
    yl = yl(find(yl==floor(yl)));
    set(ax, 'XTick', yl);
end

if ismember('y', directions)
    set(ax, 'YTickMode', 'auto');
    yl = get(ax, 'YTick');
    yl = yl(find(yl==floor(yl)));
    set(ax, 'YTick', yl);
end
