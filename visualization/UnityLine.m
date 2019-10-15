function [  ] = UnityLine( varargin )
%Makes a unity line in a plot
%Options
%   'linetype'
%   'linecolor'
%%
p = inputParser;
addParameter(p,'linetype',':')
addParameter(p,'linecolor','k')
parse(p,varargin{:})
linetype = p.Results.linetype;
linecolor = p.Results.linecolor;

%%
xrange = get(gca,'xlim');
yrange = get(gca,'ylim');

plot([min([xrange yrange]) max([xrange yrange])],...
    [min([xrange yrange]) max([xrange yrange])],linetype,'color',linecolor)

end

