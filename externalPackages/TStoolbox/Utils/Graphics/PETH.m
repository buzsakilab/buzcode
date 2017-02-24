function [fh,dArea] = PETH(S,center,TStart,TEnd, varargin)


% USAGE
%     
%  	H = PETH(S,center,TStart,TEnd, 'BinSize',binsize)
%  	S : ts
%  	center : ts of reference times
%  	TSart/TEnd : bounds of PETH
%  	optional:
%  		'BinSize' : size of bins in 10^-4s
%  Adrien Peyrache 2007

font_name = 'Arial';
    font_size = 10;
    font_weight = 'bold';
    line_width = 2;
    
opt_varargin = varargin;

defined_options  = dictArray({  { 'BinSize', {10, {'numeric'}}},
				{ 'Position', {[498   588   333   152], {'numeric'}}},
				{ 'Color', {[0 0 0], {'numeric'}}},
				{ 'ShowFigure',{1,{'numeric'}}}
                                });
getOpt;

is = intervalSet(Range(center)+TStart, Range(center)+TEnd);
sweeps = intervalSplit(S, is, 'OffsetStart', TStart);
ss = oneSeries(sweeps);
sq = intervalRate2(ss, regular_interval(TStart, TEnd, BinSize));
d = Data(sq);
sq = tsd(Range(sq),d/length(sweeps));

dArea =  Data(sq)/length(sweeps);

if ShowFigure
	
	area(Range(sq, 'ms'), Data(sq)/length(sweeps), 'FaceColor', Color);
	set(gca, 'FontName', font_name);
	set(gca, 'FontWeight', font_weight);
	set(gca, 'FontSize', font_size);
	set(gca, 'LineWidth', line_width);
	set(gca, 'XLim', [TStart TEnd]/10);
	if max(dArea) > 0
	set(gca, 'YLim', [0 max(dArea) * 1.2]);
	end
	yl = get(gca, 'YTick');
	yl = yl(find(yl==floor(yl)));
	set(gca, 'YTick', yl);
	%  set(gcf,'Position',Position)
	fh = gcf;
else
	fh = 0;
end