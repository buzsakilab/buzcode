function axh = PlotISet(IS, varargin)

opt_varargin = varargin ;

defined_options = dictArray({ {'MarkerSize', {10, {'numeric'} } }
                                              {'Marker', {'k.', {'char'} } }
                                              {'MarkerMultiColor', {0, {'numeric'} } }
                                              {'Offset', {0, {'numeric'} } }
                                              {'Spread', {1, {'numeric'} } } 
                                              {'Layers', {1, {'numeric'} } } 
                                              {'LineColor', {'k', {'char'} } }
                                              {'LineMultiColor', {0, {'numeric'} } }
                                              {'FigureHandle', {[], {'numeric'} } } 
                                              {'AxHandle', {[], {'numeric'} } }
                                              {'TUnits', {'s', {'char', 'TimeUnits'} } } 
                                            });
                                        
getOpt;

t1 = Start(IS, TUnits);
t2 = End(IS, TUnits);
nt  = length(t1);

if isempty(FigureHandle)
    FigureHandle = figure;
else
    figure(FigureHandle);
end


if ~isempty(AxHandle)
    axes(AxHandle);
end

tAll = [t1 t2 (NaN * zeros(size(t1)))];

if Layers <= 1
    y = Offset * ones(size(t1));
else
    yS = linspace(Offset-Spread, Offset+Spread, Layers);
    n = ceil(nt/Layers);
    y = repmat(yS, 1, n);
    y = y(1:nt);
    y = y(:);
end

y   = [ y y (NaN * zeros(size(t1)))];

myColors = 'kbcgyrm';
% plot marks
if MarkerMultiColor
    for i = 1:length(myColors)
        to = tAll(i:(length(myColors)):end,:);
        sz = size(to,1);
        to = reshape(to', 3*sz, 1);
        yo = y(i:(length(myColors)):end,:);
        sz = size(yo,1);
        yo = reshape(yo', 3*sz, 1);
        plot(to, yo, [myColors(i) '.'], 'MarkerSize', MarkerSize)
        hold on 
    end
else
    to = reshape(tAll', 3 * nt, 1);
    yo = reshape(y', 3 *nt, 1);
    plot(to, yo, Marker, 'MarkerSize', MarkerSize); 
end
                                            
% plot lines 

if LineMultiColor
    for i = 1:length(myColors)
        to = tAll(i:(length(myColors)):end,:);
        sz = size(to,1);
        to = reshape(to', 3*sz, 1);
        yo = y(i:(length(myColors)):end,:);
        sz = size(yo,1);
        yo = reshape(yo', 3*sz, 1);
        plot(to, yo, [myColors(i) '-']);
        hold on 
    end
    
else
    to = reshape(tAll', 3 * nt, 1);
    yo = reshape(y', 3 *nt, 1);
    plot(to, yo, '-', 'Color', LineColor); 
end

axh = gca;

