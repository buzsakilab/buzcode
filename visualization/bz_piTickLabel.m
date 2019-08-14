function [  ] = bz_piTickLabel(axis,discretization)
%bz_piTickLabel(axis) changes tick labels to 0 pi, 2pi, etc.
%
%%
if ~exist('discretization','var')
    discretization = 1;
end

switch axis
    case 'x'
        limits = xlim(gca)./pi;
    case 'y'
        limits = ylim(gca)./pi;
end
whichtick = [axis,'tick'];
whichticklabel = [axis,'ticklabel'];

ticks = round(limits(1):discretization:limits(2),1);
for ll = 1:length(ticks)
    labels{ll} = [num2str(ticks(ll)),'\pi'];
end
   
set(gca,whichtick,ticks.*pi)
set(gca,whichticklabel,labels)


end

