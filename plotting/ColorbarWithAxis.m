function [ ] = ColorbarWithAxis(cbounds,label,varargin)
%ColorbarWithAxis(cbounds,label) adds a colorbar to a imageplot and labels
%axis and bounds with > and <.
%
%TO DO
%   change inputs to 'parms' with inputParser
%   add 3std functionality
%
%DLevenstein 2016
%% inputParse for input options

labelloc = 'side';
scale = 'lin';

barlocation = 'eastoutside';

%%
cb = colorbar(barlocation);

if isstring(cbounds); switch cbounds
        case '3std'
            display('ColorbarWithAxis does not yet have this functionality, don"t you wish it did?...')
        case 'datarange'
            display('ColorbarWithAxis does not yet have this functionality, don"t you wish it did?...')
end; end
            

caxis(cbounds);
switch labelloc
    case 'side'
        ylabel(cb, label);
    case 'top'
        title(cb, label);
end

cb.Ticks = [cbounds(1) mean(cbounds) cbounds(2)];

switch scale
    case 'log'
        cb.TickLabels = {['< ',num2str(2.^(cbounds(1)))],2.^(mean(cbounds)),...
            ['> ',num2str(2.^(cbounds(2)))]};
    case 'lin'
        cb.TickLabels = {['< ',num2str(cbounds(1))],mean(cbounds),...
            ['> ',num2str(cbounds(2))]};
end

end

