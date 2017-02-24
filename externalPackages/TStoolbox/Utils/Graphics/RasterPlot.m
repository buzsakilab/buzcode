function RasterPlot(S, varargin)
% RasterPlot(S, height, bar, tstart, tend)
% 
% Plot spike trains as rasters.
% INPUT:
% S: tsdArray of ts objects, spike trains
% OPTIONS: 
% Height: the height allocated to each cell (includes spacing), default 1
% BarFraction: the fraction of height occupied by a cell, default 0.8
% TStart, TEnd: the beginning of end of time interval to be plot, default
% 'Markers': times denoting a special event for each trial which that will
% be marked by a symbol
% 'MarkerTypes': symobels to denote markers, (first character color
% (optional), second character symbol. Defaults to '*'
% the whole thing
  
% batta 2001 under construction
  

opt_varargin = varargin;

defined_options = dictArray({ {'Height', {1, {'numeric'} } }, 
                                              {'BarFraction', {0.8, {'numeric'} } },
                                              {'BarColor', {{'k'}, {'cell'} } },
                                               {'LineWidth', {1, {'numeric'} } },
                                              {'TStart', {NaN, {'numeric'} } },
                                              {'TEnd', {1e100, {'numeric'} } } 
                                              { 'FigureHandle', {[], {'numeric'} } } 
                                              { 'AxHandle', {[], {'numeric'} } } 
                                              {'Markers', {{}, {'cell'} } }
                                              {'MarkerTypes', { {}, {'cell' } } }
                                              {'MarkerSize', { [], {'numeric'} } }
					      {'MarkerEdgeColor', {{'r'}, {'cell'} } }
					      {'MarkerFaceColor', {{''}, {'cell'} } }
                                              {'TimeUnits', {{'ms'}, {'cell'} } }
                                              {'YOffset', {0, {'numeric'} } }
                                              {'XOffset', {0, {'numeric'} } }
                                           });

getOpt;


if isempty(FigureHandle)
    FigureHandle = gcf;
else
    figure(FigureHandle);
end

% 
% if ~isempty(AxHandle)
%     
%     axes(AxHandle);
% end

if isfinite(TStart )
  for i = 1:length(S)
    Si{i} = Restrict(S{i}, TStart, TEnd);
  end
else 
  Si = S;
end


  
for i = 1:length(Si)
  sp = Range(Si{i}, TimeUnits{1});
  sx = [sp sp repmat(NaN, length(sp), 1)];
  sy = repmat([(i*Height) (i*Height + Height *BarFraction) NaN], length(sp), 1);
  sx = reshape(sx', 1, length(sp)*3);
  sy = reshape(sy', 1, length(sp)*3);
  
  line(sx-XOffset, sy+YOffset, 'Color', BarColor{1}, 'LineWidth', LineWidth);
%    set(gca, 'ylim', [1 length(Si)+1]);
  
  hold on

end



if ~isempty(Markers)
    if isempty(MarkerTypes)
        for iM  = 1:length(Markers)
            MarkerTypes{iM} = '*';
        end
    end
    
    for iM = 1:length(Markers)
        M = Markers{iM};
        mt = MarkerTypes{iM};
        if length(mt) == 2
            mr = mt(2);
            mc = mt(1);
	    mf = mt(1);
        elseif length(mt)==1
            mr = mt;
            mc= MarkerEdgeColor{1};
            mf= MarkerFaceColor{1};
        else error('MarkerTypes elements should have one or two characters')
        end
        
        i = 1:length(Si);
	if length(mf)
		l = line(M, (i+1/2)*Height, 'LineStyle', 'none', 'Marker', mr, 'MarkerEdgeColor', mc,'MarkerFaceColor',mf);
	else
		l = line(M, (i+1/2)*Height, 'LineStyle', 'none', 'Marker', mr, 'MarkerEdgeColor', mc);
	end
        if ~isempty(MarkerSize)
            set(l, 'MarkerSize', 20);
        end
    end
end


yl = get(gca, 'YTick');
yl = yl(find(yl==floor(yl)));
set(gca, 'YTick', yl);