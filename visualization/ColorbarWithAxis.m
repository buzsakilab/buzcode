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
parms = inputParser;
addParameter(parms,'location','eastoutside');
addParameter(parms,'labelloc','side');
addParameter(parms,'scale','lin');
addParameter(parms,'inclusive',{'',''});


parse(parms,varargin{:})
barlocation = parms.Results.location;
labelloc = parms.Results.labelloc;
scale = parms.Results.scale;
inclusive = parms.Results.inclusive;

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
        cb.TickLabels = {[inclusive{1},num2str(2.^(cbounds(1)))],2.^(mean(cbounds)),...
            [inclusive{2},num2str(2.^(cbounds(2)))]};
    case 'lin'
        cb.TickLabels = {[inclusive{1},num2str(cbounds(1))],mean(cbounds),...
            [inclusive{2},num2str(cbounds(2))]};
end

end

