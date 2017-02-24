function [fh, rasterAx] = Raster(S, TStart, TEnd, varargin)



font_name = 'Arial';
    font_size = 10;
    font_weight = 'bold';
    line_width = 2;
    
opt_varargin = varargin;

defined_options  = dictArray({  { 'BinSize', {10, {'numeric'}}},
                                {'LineWidth', {1, {'numeric'} } },
                                {'ShowXTick', {1, {'numeric'} } },
                                {'ShowXTickLabel', {1, {'numeric'} } },
                                {'Ax', {0, {'numeric'}}}  });
getOpt;

is = intervalSet(TStart, TEnd);
sweeps = Restrict(S,is);
    fh = gcf;

if ~Ax
    rasterAx = axes('position', [0.1 0.05 0.8 0.8]);
    fh = gcf;
    axes(rasterAx);
else
    rasterAx = Ax;
end

set(gca, 'FontName', font_name);
set(gca, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', line_width);
set(gca, 'XLim', [TStart TEnd]/10);
set(gca, 'XTick', [TStart TEnd]/10);
if ~ShowXTick
	set(gca,'XTick',[]);
end
if ~ShowXTickLabel
	set(gca,'XTickLabel',{});
end

RasterPlot(sweeps, 'AxHandle', rasterAx, ...
    'FigureHandle', fh, ...
    'TStart', TStart, ...
    'TEnd', TEnd, ...
    'LineWidth', LineWidth);

fh = gcf;
