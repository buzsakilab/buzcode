classdef (CaseInsensitiveProperties = true) functionalSpreadPlot < ...
        iosr.statistics.functionalPlot
%FUNCTIONALSPREADPLOT Draw a functional plot showing data spread
% 
%   Use this class to plot the spread of data against a continuous
%   variable. The class calculates point-wise statistics for each set of
%   observations on the x-axis. For a plot of spread where each _function_
%   constitutes a single observation, use IOSR.STATISTICS.FUNCTIONALBOXPLOT
%   instead.
% 
%   Measures of spread include parametric and non-parametric confidence
%   intervals, interquartile range, standard deviation, or arbitrary
%   percentile ranges. The spread is indicated by a shaded region, with a
%   central line indicating the primary metric (mean or median). Additional
%   percentiles, outliers, and the limits of the data may also be plotted.
% 
%   FUNCTIONALSPREADPLOT operates on the array Y, which is an N-by-M-by-P
%   array, where N observations for each X value are stored in the columns
%   of Y, there are M X-values, and there are P functional plots.
%   
%   IOSR.STATISTICS.FUNCTIONALSPREADPLOT properties:
%       addPrctiles             - Plot lines for additional percentiles.
%                                 The property should be a vector of
%                                 percentiles; each percentile will be
%                                 plotted for each functional plot. The
%                                 property is empty by default.
%       addPrctilesLineColor    - Specify the line colors for the
%                                 additional percentile lines. The property
%                                 may be specified using any valid MATLAB
%                                 color specifier. Multiple colors (for
%                                 multiple function plots) may be specified
%                                 by passing colors in a length-P cell
%                                 vector, or by passing a color map
%                                 function handle. The default is 'auto',
%                                 which uses Matlab's 'line' color map.
%       addPrctilesLineStyle    - Specify the line styles for the
%                                 additional percentile lines. The property
%                                 may be specified as a char array (e.g.
%                                 '-'), or as a P-by-Q cell array of
%                                 strings, where Q is the number of
%                                 additional percentiles. The default is
%                                 '--'.
%       addPrctilesLineWidth    - Specify the line widths for the
%                                 additional percentile lines. The property
%                                 may be specified as a numeric array (e.g.
%                                 '-'), or as a P-by-Q cell array, where Q
%                                 is the number of additional percentiles.
%                                 The default is 1.
%       innerMode               - Specify the mode for plotting the inner
%                                 shaded area of the plot. The options are:
%                                   - 'quartile' : the 25th to 75th
%                                                  percentile
%                                                  (interquartile range)
%                                   - 'std'      : the standard deviation
%                                   - 'notch'    : the non-parametric
%                                                  confidence interval (see
%                                                  iosr.statistics.boxPlot)
%                                   - 'ci90'     : the 90% confidence
%                                                  interval
%                                   - 'ci95'     : the 95% confidence
%                                                  interval
%                                   - 'ci99'     : the 99% confidence
%                                                  interval
%                                 The default is 'quartile'.
%       limit                   - Specifies the mode indicating the limits
%                                 that define outliers. See 
%                                 iosr.statistics.boxPlot. The default is
%                                 '1.5IQR'.
%       mainMode                - Specifies the statistic used for the
%                                 central line. The options are 'median' or
%                                 'mean'. The default is 'median'.
%       method                  - Specifies the method used to calculate
%                                 the quantiles. See
%                                 iosr.statistics.boxPlot. The default is
%                                 'R-8'.
%       outerMode               - Specifies what the outer lines plot. The
%                                 options are 'limit' (plots the limits of
%                                 the data as determined by the 'limit'
%                                 option) or 'none'.
%       outlierSize             - Specifies the size, in square points, of
%                                 the outlier markers. The default is 36.
%       outlierEdgeColor        - Specifies the color of the outlier marker
%                                 edges. See 'addPrctilesLineColor' for
%                                 details of how to specify colors. The
%                                 default is 'auto'.
%       outlierFaceColor        - Specifies the color of the outlier marker
%                                 faces. See 'addPrctilesLineColor' for
%                                 details of how to specify colors. The
%                                 default is 'none'.
%       outlierMarker           - Specifies the markers used for the
%                                 outliers. The property may be specified
%                                 as a char array (e.g. '+'), or as a cell
%                                 array of strings (if passing markers for
%                                 each functional plot). The default is
%                                 '+'.
% 
%   Additional properties are inherited from
%   IOSR.STATISTICS.FUNCTIONALPLOT.
% 
%   These properties can be referenced using dot notation - e.g. H.LIMIT
%   where H is an instance of the FUNCTIONALSPREADPLOT object - or using
%   the SET and GET methods - e.g. GET(H,'LIMIT'). Both methods are
%   case-insensitive.
% 
%   Note that some handles will be empty unless the associated option is
%   turned on.
% 
%   Read-only properties:
%       x                   - The x data.
%       y                   - The y data.
%       weights             - Array giving the weights for the data in y.
%                             The array must be the same size as y.
%       statistics          - Structure containing the statistics used
%                             for the box plot. With the exception of
%                             'outliers', 'outliers_IX', 'addPrctiles', and
%                             'percentile' noted below, each field contains
%                             a 1-by-M-by-P numeric array of values
%                             (identical to that returned by
%                             IOSR.STATISTICS.QUANTILE). The fields are:
%                                 'addPrctiles'   : an A-by-M-by-P... array
%                                                   containing the
%                                                   additional percentile
%                                                   values specified by the
%                                                   'addPrctiles' property,
%                                                   where A is the number
%                                                   of percentiles in the
%                                                   property
%                                 'percentile'    : the percentile limits
%                                                   specified by the
%                                                   'percentile' property
%                                 'median'        : the median values
%                                 'N'             : the sample size
%                                 'Q1'            : the 25th percentiles
%                                 'Q3'            : the 75th percentiles
%                                 'inner_l'       : the inner shaded area's
%                                                   lower limit
%                                 'inner_u'       : the inner shaded area's
%                                                   upper limit
%                                 'IQR'           : the inter-quartile
%                                                   ranges
%                                 'main'          : the values of the
%                                                   central line
%                                 'mean'          : the mean values
%                                 'min'           : the minimum values
%                                                   (excl. outliers)
%                                 'max'           : the maximum values
%                                                   (excl. outliers)
%                                 'notch_u'       : the upper notch values
%                                 'notch_l'       : the lower notch values
%                                 'outer_l'       : the outer line's lower
%                                                   limit
%                                 'outer_u'       : the outer line's upper
%                                                   limit
%                                 'outliers'      : a 1-by-M-by-P cell
%                                                   array of outlier values
%                                 'outliers_IX'   : a logical array, the
%                                                   same size as Y, with
%                                                   true values indicating
%                                                   outliers
%                                 'std'          : the standard deviations
% 
%   IOSR.STATISTICS.FUNCTIONALSPREADPLOT methods:
%       functionalSpreadPlot    - Create the plot.
% 
%   See also IOSR.STATISTICS.BOXPLOT, IOSR.STATISTICS.QUANTILE,
%   IOSR.STATISTICS.FUNCTIONALPLOT, IOSR.STATISTICS.FUNCTIONALBOXPLOT.

    properties (AbortSet)
        addPrctiles = []                % Additional percentiles to plot.
        addPrctilesLineColor = 'auto'   % Colors for the additional percentile lines.
        addPrctilesLineStyle = '--'     % Styles for the additional percentile lines.
        addPrctilesLineWidth = 0.5      % Widths for the additional percentile lines.
        innerMode = 'quartile'          % Mode of the central shaded area.
        limit = '1.5IQR'                % Mode indicating the limits that define outliers.
        mainMode = 'median'             % The statistic used for the central line.
        method = 'R-8'                  % The method used to calculate the quantiles.
        outerMode = 'limit'             % Mode of the outer lines.
        outlierSize = 36                % Size of the outlier markers.
        outlierEdgeColor = 'auto'       % Color of the outlier marker edges.
        outlierFaceColor = 'none'       % Color of the outlier marker faces.
        outlierMarker = '+'             % The outlier marker.
    end
    
    properties (Access = private)
        data
    end
    
    methods
        
        % constructor
        function obj = functionalSpreadPlot(varargin)
        % FUNCTIONALSPREADPLOT Draw a functional spread plot.
        %
        %   IOSR.STATISTICS.FUNCTIONALSPREADPLOT(Y) produces a functional
        %   spread plot of the data in Y. Y should be an N-by-M or
        %   N-by-M-by-P array, where N observations for each X value are
        %   stored in the columns of Y, there are M X-values, and there are
        %   P functional plots. The plot shows the median as the main line,
        %   the interquartile range as the shaded region, lines showing the
        %   limits of the data, and markers showing any outliers. See
        %   IOSR.STATISTICS.BOXPLOT for an explanation of how limits and
        %   outliers are estimated. This plot is traditionally known as a
        %   functional box plot.
        %
        %   Tabular data can be arranged into the appropriate format using
        %   the IOSR.STATISTICS.TAB2BOX function.
        %
        %   IOSR.STATISTICS.FUNCTIONALSPREADPLOT(X,Y) specifies the x-axis
        %   values for the plot. X should be an M-length vector. The
        %   default is 1:M.
        %
        %   IOSR.STATISTICS.FUNCTIONALSPREADPLOT(...,'PARAMETER',VALUE)
        %   allows the plotting options to be specified when the plot is
        %   constructed.
        %
        %   IOSR.STATISTICS.FUNCTIONALSPREADPLOT(AX,...) creates the plot
        %   in the axes specified by AX.
        %   
        %   Examples
        %
        %     Example 1: Draw a functional box plot
        %       % generate random data
        %       y = cat(3, randn(100,21), 0.25+randn(100,21));
        %       x = 0:20;
        %    
        %       % Draw a functional box plot for the first function
        %       figure
        %       iosr.statistics.functionalSpreadPlot(x, y(:,:,1));
        %       title('Functional box plot.')
        %       axis tight
        %       box on
        % 
        %     Example 2: Plot means and confidence intervals
        %       figure
        %       iosr.statistics.functionalSpreadPlot(x, y,...
        %           'spreadAlpha',0.5,...
        %           'outerMode', 'none', ...
        %           'innermode', 'ci95', ...
        %           'mainmode', 'mean', ...
        %           'showOutliers', false);
        %       title('Mean and 95% confidence intervals')
        %       axis tight
        %       box on
        
            start = obj.getXY(varargin{:});
            
            % check input is valid size
            assert(ndims(obj.y) <= 3 && ndims(obj.y) >= 2, ...
                'iosr:functionalSpreadPlot:invalidY', ...
                'Y must be a two- or three-dimensional array.');
            
            % check weights
            obj.checkWeights();
            
            % set the properties of the plot
            obj.setProperties(start,nargin,varargin);
            
            % remove NaN columns
            obj.removeNaN();
            
            % calculate the statistics used in the plot
            obj.calculateStats();
            
            %% draw
            
            % set handles
            obj.parseAxesHandle(varargin{:});
            
            % draw the box plot
            obj.draw();
            
        end
        
        %% Accessors / Mutators
        
        function set.mainMode(obj, val)
            obj.mainMode = val;
            obj.calculateStats();
            obj.draw();
        end
        
        function set.innerMode(obj, val)
            obj.innerMode = val;
            obj.calculateStats();
            obj.draw();
        end
        
        function set.limit(obj, val)
            obj.limit = val;
            obj.calculateStats();
            obj.draw();
        end
        
        function set.method(obj, val)
            obj.method = val;
            obj.calculateStats();
            obj.draw();
        end
        
        function set.outerMode(obj, val)
            obj.outerMode = val;
            obj.calculateStats();
            obj.draw();
        end
        
        % additional percentiles
        
        function set.addPrctiles(obj, val)
            if ~isempty(val)
                assert(all(val >= 0) && all(val <= 1000), 'iosr:functionalSpreadPlot:invalidPrctiles', ...
                    'Additional percentiles must be in the range [0, 100].')
                obj.addPrctiles = val;
            end
            obj.calculateStats();
            obj.draw();
        end
        
        function val = get.addPrctilesLineColor(obj)
            val = obj.parseColor(obj.addPrctilesLineColor);
        end
        
        function set.addPrctilesLineColor(obj,val)
            obj.addPrctilesLineColor = val;
            obj.draw();
        end
        
        function val = get.addPrctilesLineStyle(obj)
            val = obj.parseProps(obj.addPrctilesLineStyle, true);
        end
        
        function set.addPrctilesLineStyle(obj,val)
            obj.addPrctilesLineStyle = val;
            obj.draw();
        end
        
        function val = get.addPrctilesLineWidth(obj)
            val = obj.parseProps(obj.addPrctilesLineWidth, true);
        end
        
        function set.addPrctilesLineWidth(obj,val)
            obj.addPrctilesLineWidth = val;
            obj.draw();
        end
        
        % outlier settings

        function set.outlierSize(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:functionalSpreadPlot:invalidOutlierSize', ...
                '''OUTLIERSIZE'' must be a numeric scalar')
            obj.outlierSize = val;
            obj.draw();
        end
        
        function val = get.outlierEdgeColor(obj)
            val = obj.parseColor(obj.outlierEdgeColor);
        end
        
        function set.outlierEdgeColor(obj,val)
            obj.outlierEdgeColor = val;
            obj.draw();
        end
        
        function val = get.outlierFaceColor(obj)
            val = obj.parseColor(obj.outlierFaceColor);
        end
        
        function set.outlierFaceColor(obj,val)
            obj.outlierFaceColor = val;
            obj.draw();
        end
        
        function val = get.outlierMarker(obj)
            val = obj.parseProps(obj.outlierMarker, false);
        end
        
        function set.outlierMarker(obj,val)
            obj.outlierMarker = val;
            obj.draw();
        end
        
    end
    
    methods (Access = protected)
        
        % calculate the statistics used in the plot
        function calculateStats(obj)
            
            % Most statistics are calculated in the base class
            calculateStats@iosr.statistics.statsPlot(obj);

            % Calculate specific stats
            outsize = size(obj.statistics.median);
            obj.statistics.addPrctiles = zeros([length(obj.addPrctiles) outsize(2:end)]);
            obj.statistics.main = zeros(outsize);
            obj.statistics.inner_u = zeros(outsize);
            obj.statistics.inner_l = zeros(outsize);
            obj.statistics.outer_u = zeros(outsize);
            obj.statistics.outer_l = zeros(outsize);
            subidx = cell(1,length(obj.outDims));
            for n = 1:prod(obj.outDims)
                [subidx{:}] = ind2sub(obj.outDims, n);
                subidxAll = subidx;
                subidxAll{1} = ':';

                if isnumeric(obj.innerMode)
                    assert(numel(obj.innerMode) == 2, 'iosr:functionalSpreadPlot:invalidInnerMode', ...
                        'If numeric, INNERMODE must be a two-element vector.')
                    obj.statistics.inner_u(subidx{:}) = iosr.statistics.quantile(obj.y(subidxAll{:}), ...
                        max(obj.innerMode)/100, [], obj.method, obj.weights(subidxAll{:}));
                    obj.statistics.inner_l(subidx{:}) = iosr.statistics.quantile(obj.y(subidxAll{:}), ...
                        min(obj.innerMode)/100, [], obj.method, obj.weights(subidxAll{:}));
                else
                    assert(ischar(obj.innerMode), 'iosr:functionalSpreadPlot:invalidInnerMode', ...
                        'INNERMODE must be numeric or a char array.')
                    switch obj.innerMode
                        case 'quartile'
                            obj.statistics.inner_u(subidx{:}) = obj.statistics.Q3(subidx{:});
                            obj.statistics.inner_l(subidx{:}) = obj.statistics.Q1(subidx{:});
                        case 'notch'
                            obj.statistics.inner_u(subidx{:}) = obj.statistics.notch_u(subidx{:});
                            obj.statistics.inner_l(subidx{:}) = obj.statistics.notch_l(subidx{:});
                        case 'std'
                            std = obj.statistics.std(subidx{:})/2;
                            obj.statistics.inner_u(subidx{:}) = obj.statistics.mean(subidx{:}) + std;
                            obj.statistics.inner_l(subidx{:}) = obj.statistics.mean(subidx{:}) - std;
                        case 'ci90'
                            cil = ( (1.645 * obj.statistics.std(subidx{:})) / sqrt(obj.statistics.N(subidx{:})) );
                            obj.statistics.inner_u(subidx{:}) = obj.statistics.mean(subidx{:}) + cil;
                            obj.statistics.inner_l(subidx{:}) = obj.statistics.mean(subidx{:}) - cil;
                        case 'ci95'
                            cil = ( (1.96 * obj.statistics.std(subidx{:})) / sqrt(obj.statistics.N(subidx{:})) );
                            obj.statistics.inner_u(subidx{:}) = obj.statistics.mean(subidx{:}) + cil;
                            obj.statistics.inner_l(subidx{:}) = obj.statistics.mean(subidx{:}) - cil;
                        case 'ci99'
                            cil = ( (2.576 * obj.statistics.std(subidx{:})) / sqrt(obj.statistics.N(subidx{:})) );
                            obj.statistics.inner_u(subidx{:}) = obj.statistics.mean(subidx{:}) + cil;
                            obj.statistics.inner_l(subidx{:}) = obj.statistics.mean(subidx{:}) - cil;
                        otherwise
                            error('iosr:functionalSpreadPlot:unknownInnerMode','Unkonwn INNERMODE.')
                    end
                end

                if isnumeric(obj.outerMode)
                    assert(numel(obj.outerMode) == 2, 'iosr:functionalSpreadPlot:invalidOuterMode', ...
                        'If numeric, OUTERMODE must be a two-element vector.')
                    obj.statistics.outer_u(subidx{:}) = iosr.statistics.quantile(obj.y(subidxAll{:}), ...
                        max(obj.outerMode)/100, [], obj.method, obj.weights(subidxAll{:}));
                    obj.statistics.outer_l(subidx{:}) = iosr.statistics.quantile(obj.y(subidxAll{:}), ...
                        min(obj.outerMode)/100, [], obj.method, obj.weights(subidxAll{:}));
                else
                    assert(ischar(obj.outerMode), 'iosr:functionalSpreadPlot:invalidOuterMode', ...
                        'OUTERMODE must be numeric or a char array.')
                    switch obj.outerMode
                        case 'limit'
                            obj.statistics.outer_u(subidx{:}) = obj.statistics.max(subidx{:});
                            obj.statistics.outer_l(subidx{:}) = obj.statistics.min(subidx{:});
                        case 'none'
                            obj.statistics.outer_u(subidx{:}) = NaN;
                            obj.statistics.outer_l(subidx{:}) = NaN;
                        otherwise
                            error('iosr:functionalSpreadPlot:unknownOuterMode','Unknown OUTERMODE.')
                    end
                end

                switch obj.mainMode
                    case 'mean'
                        obj.statistics.main(subidx{:}) = obj.statistics.mean(subidx{:});
                    case 'median'
                        obj.statistics.main(subidx{:}) = obj.statistics.median(subidx{:});
                    otherwise
                            error('iosr:functionalSpreadPlot:unknownMainMode','Unknown MAINMODE.')
                end

                if ~isempty(obj.addPrctiles)
                    obj.statistics.addPrctiles(subidxAll{:}) = ...
                        iosr.statistics.quantile(obj.y(subidxAll{:}),obj.addPrctiles(:)./100,[],obj.method,obj.weights(subidxAll{:}));
                end

                temp = obj.y(subidxAll{:});
                obj.statistics.mean(subidx{:}) = mean(temp(~isnan(temp)));
            end
            
        end
        
        % Draw the plot
        function draw(obj)
            
            if isfield(obj.handles,'axes')

                axes(obj.handles.axes);
                cla;
                hold on;

                if length(obj.outDims) > 2
                    nlines = obj.outDims(3);
                else
                    nlines = 1;
                end
                
                draw@iosr.statistics.functionalPlot(obj);

                for n = nlines:-1:1

                    % draw additional percentiles
                    if ~isempty(obj.addPrctiles)
                        for p = 1:numel(obj.addPrctiles)
                            obj.handles.addPrctiles(p,n) = line(obj.x, squeeze(obj.statistics.addPrctiles(p,:,n)), ...
                                'color', obj.addPrctilesLineColor{n}, ...
                                'linewidth', obj.addPrctilesLineWidth{n,p}, ...
                                'linestyle', obj.addPrctilesLineStyle{n,p});
                        end
                    end

                    % draw outliers
                    if (obj.showOutliers)
                        xOutliers = repmat(obj.x,size(obj.y, 1), 1);
                        xOutliers = xOutliers(obj.statistics.outliers_IX(:,:,n));
                        yOutliers = obj.y(:,:,n);
                        yOutliers = yOutliers(obj.statistics.outliers_IX(:,:,n));
                        obj.handles.outliers(n) = scatter(xOutliers, yOutliers, ...
                            obj.outliersize, obj.outlierMarker{n}, ...
                            'MarkerEdgeColor', obj.outlierEdgeColor{n}, ...
                            'MarkerFaceColor', obj.outlierFaceColor{n}, ...
                            'LineWidth', obj.outlierLineWidth{n});
                    end

                end

            end
            
        end
        
    end
    
end
