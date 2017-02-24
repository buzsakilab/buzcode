function [fh, rasterAx, histAx,matVal] = ImagePETH2(S, center, TStart, TEnd, varargin)

% Display a PETH-like plot for continuous variable such as EEG
% Uses RasterImagePlot2
% Adrien Peyrache 2007, adapted from Battglia


font_name = 'Arial';
    font_size = 10;
    font_weight = 'bold';
    line_width = 2;
    figure(2),clf
caxis([-73 -55]) 
cb = colorbar('Position',[0.8 0.1 0.04 0.3]);
set(cb,'YTick',[0 63])
set(cb,'YTickLabel',[-73 -55])

opt_varargin = varargin;

defined_options  = dictArray({ { 'RasterFraction', {0.7, {'numeric'}} }
                               { 'BinSize', {10, {'numeric'}}},
                                {'LineWidth', {1, {'numeric'} } },
                                {'Markers', { {}, {'cell'}} } ,
                                 {'MarkerTypes', { {}, {'cell'}}}, 
                                {'MarkerSize', { [], {'numeric'} } }
                                });
getOpt;

is = intervalSet(Range(center)+TStart, Range(center)+TEnd);

sweeps = intervalSplit(S, is, 'OffsetStart', TStart);

for iM = 1:length(Markers)
    Markers{iM} = (Range(Markers{iM}) - Range(center))/10000; 
end

rf = RasterFraction * 0.8;
rasterAx = axes('position', [0.1 0.05 0.8 (rf+0.05)]);
histAx = axes('position', [0.1 (rf+0.15) 0.8 (0.75-rf)]);

fh = gcf;
axes(rasterAx);


set(gca, 'FontName', font_name);
set(gca, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', line_width);
set(gca, 'XLim', [TStart TEnd]/10000);
matVal = RasterImagePlot(sweeps, 'AxHandle', rasterAx, ...
	'FigureHandle', fh, ...
	'TStart', TStart, ...
	'TEnd', TEnd, ...
	'LineWidth', LineWidth, ...
	'Markers', Markers, ...
	'MarkerTypes', MarkerTypes, ...
	'MarkerSize', MarkerSize);

set(gca, 'Box', 'on');
axes(histAx);


%  keyboard

dArea =  mean(Data(matVal)');

%afrea(Range(matVal, 's'), dArea, 'FaceColor', 'k');
plot(Range(matVal, 's'), dArea);

set(gca, 'FontName', font_name);
set(gca, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', line_width);
set(gca, 'XLim', [TStart TEnd]/10000);
% 
% yl = [min(0,1.2*min(dArea)) max(dArea) * 1.2];
% 
% if max(dArea) > 0
%     set(gca, 'YLim', yl);
% end
% 
% yM = max(yl);%floor(100*max(yl))/100;
% ym = min(yl);%floor(100*min(yl))/100;
% yl = [ym:(yM-ym)/3:yM];
% set(gca, 'YTick', yl);
fh = gcf;
