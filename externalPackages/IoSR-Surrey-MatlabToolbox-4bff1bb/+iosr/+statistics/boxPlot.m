classdef (CaseInsensitiveProperties = true) boxPlot < iosr.statistics.statsPlot
%BOXPLOT Draw a box plot
%   
%   Use this class to draw box plots. The class provides a number of
%   options not included in Matlab's boxplot function.
% 
%   For each box, the central mark is the median, the box extends
%   vertically between the 25th and 75th percentiles, the whiskers extend
%   to the most extreme data that are not considered outliers, and the 
%   outliers are plotted individually.
% 
%   BOXPLOT can draw boxes for data in Y for an arbitrary number of
%   dimensions. If Y is an N-by-P-by-G-by-I-by-J... array then G*I*J*...
%   boxes are plotted hierarchically for each column P; i.e. J boxes are
%   plotted for each index of I, and J*I boxes are plotted for each index
%   of G.
%   
%   IOSR.STATISTICS.BOXPLOT properties:
%     Plotting options:
%       addPrctiles         - Show additional percentiles using markers and
%                             labels. The property should be a vector of
%                             percentiles; each percentile will be plotted
%                             for each box. The property is empty by
%                             default.
%       addPrctilesColor    - Specify the marker color for the additional
%                             percentile markers. The property should be a
%                             cell array of strings indicating the color
%                             for each percentile; the colors will be
%                             repeated for each box. The default is black.
%       addPrctilesLabels   - Specify labels for the additional
%                             percentiles. The property should be a cell
%                             array of strings. By defualt no label
%                             is shown.
%       addPrctilesMarkers  - Specify markers for the additional
%                             percentiles. The property should be a cell
%                             array of strings indicating the shape of each
%                             percentile; the markers will be repeated for
%                             each box. The default is '*'.
%       addPrctilesSize     - Specify the marker size for the additional
%                             percentile markers. The property should be a
%                             numeric vector indicating the size for each
%                             percentile; the markers will be repeated for
%                             each box. The default is 6.
%       addPrctilesTxtSize  - Specify the font size of the additional
%                             percentile labels. The property should be a 
%                             numeric scalar. The default is 9.
%       boxAlpha            - The transparency of the boxes (0 is
%                             transparent, 1 is opaque). The default is 1.
%       boxColor            - Fill color of the box. The setting can be a
%                             color specifier (3-element vector or single
%                             character), 'auto', 'none', or a handle to a
%                             colormap function (e.g. @gray); 'auto' means
%                             that Matlab's default colors are used. The
%                             default is 'none'.
%       boxWidth            - The width of the box. The setting can be
%                             'auto' (meaning that the width is determined
%                             automatically) or a scalar (specifying the
%                             width in x-axis units). The default is
%                             'auto'.
%       groupLabelFontSize  - Specify the font size of the group labels
%                             when the 'style' option is set to
%                             'hierarchy'. The default is 9.
%       groupLabelHeight    - Specify the height of the area reserved for
%                             group labels when the 'style' option is set
%                             to 'hierarchy'. The height is determined
%                             automatically ('auto') or can be specified in
%                             scale units. The default is 'auto'.
%       groupLabels         - If the 'style' option is set to 'hierarchy'
%                             then these labels will be used to label the
%                             boxes for each x-group. The parameter should
%                             be specified as a cell vector whereby the Nth
%                             element contains a vector of length
%                             SIZE(Y,N+2). By default, the labels are
%                             empty.
%       groupWidth          - Specify the proportion of the x-axis interval
%                             across which each x-group of boxes should be
%                             spread. The default is 0.75.
%       limit               - Mode indicating the limits that define
%                             outliers. When set to '1.5IQR' or '3IQR', the
%                             min and max values are Q1-Z*IQR and Q3+Z*IQR,
%                             where Z = 1.5 or 3 respectively. When set to
%                             'none', the min and max values are the min
%                             and max of the data (in this case there will
%                             be no outliers). The property may also be
%                             specified as a two-element vector determining
%                             the limits as percentiles (e.g. [5,95]). The
%                             default is '1.5IQR'.
%       lineColor           - Color of the box outline and whiskers. The
%                             default is 'k'.
%       lineStyle           - Style of the whisker line. The default is
%                             '-'.
%       lineWidth           - Width, in points, of the box outline, whisker
%                             lines, notch line, violin, and outlier marker
%                             edges. The default is 1.
%       meanColor           - Color of the mean marker when showMean=true.
%                             The default is 'auto' (see boxColor).
%       meanMarker          - Marker used for the mean when showMean=true.
%                             The default is '+'.
%       meanSize            - Size, in point units, of the mean markers
%                             when showMean=true. The default is 6.
%       medianColor         - Color of the median line. The default is
%                             'auto' (see boxColor).
%       method              - The method used to calculate the quantiles,
%                             labelled according to
%                             http://en.wikipedia.org/wiki/Quantile. The
%                             default is 'R-8' whereas the default for
%                             Matlab is 'R-5'. See
%                             IOSR.STATISTICS.QUANTILE.
%       notch               - Logical value indicating whether the box
%                             should have a notch. The notch is centred on
%                             the median and extends to ±1.58*IQR/sqrt(N),
%                             where N is the sample size (number of non-NaN
%                             rows in Y). Generally if the notches of two
%                             boxes do not overlap, this is evidence of a
%                             statistically significant difference between
%                             the medians. The default is false.
%       notchDepth          - Depth of the notch as a proportion of half
%                             the box width. The default is 0.4.
%       notchLine           - Logical value specifying whether to draw a
%                             horizontal line in the box at the extremes of
%                             the notch. (May be specified indpendently of
%                             'notch'.) The default is false.
%       notchLineColor      - Color of the notch line when notchLine=true.
%                             The default is 'k'.
%       notchLineStyle      - Line style of the notch line when
%                             notchLine=true. The default is ':'.
%       outlierSize         - Size, in square points, of the outlier
%                             marker. The default is 36.
%       percentile          - Percentile limits of the boxes. The default
%                             is [25,75]. All other statistics, such as
%                             those that depend on the IQR (e.g. whisker
%                             extent and notch height), are unaffected by
%                             this parameter.
%       sampleFontSize      - Specify the font size of the sample size
%                             display (if sampleSize=true). The default is
%                             9.
%       sampleSize          - Specify whether to display the sample size in
%                             the top-right of each box. The default is
%                             false.
%       scaleWidth          - Logical value to specify scaling the width of
%                             each box according to the square root of the
%                             sample size. The default is false.
%       scatterAlpha        - Set the transparency of the scatter markers
%                             (0 is transparent, 1 is opaque) when
%                             showScatter=true. The default is 1.
%       scatterColor        - Scatter marker color for scatter plots of
%                             underlying data (if showScatter=true). The
%                             default is [.5 .5 .5].
%       scatterLayer        - Set the layer of scatter plots with respect
%                             to the boxes. Can be set to 'top' or
%                             'bottom'. The default is 'top'.
%       scatterMarker       - Marker used for scatter plots of underlying
%                             data (if showScatter=true). The default is
%                             'x'.
%       scatterSize         - Size, in square points, of the scatter
%                             markers (if showScatter=true). The default is
%                             36.
%       showLegend          - Display a legend of the data. The labels can
%                             be set using the 'groupLabels' option. Note
%                             that the legend uses the box color, median
%                             line, or violin fill color to distinguish
%                             legend entries. If these properties do not
%                             differ between boxes then an error will be
%                             returned.
%       showMean            - Logical value determines whether to display
%                             the mean of the data for each box. The
%                             default is false.
%       showOutliers        - Logical value determines whether to display
%                             outliers. The default is true.
%       showScatter         - If set to true, a scatter plot of the
%                             underlying data for each box will be
%                             overlayed on the plot. The data will have a
%                             random x-axis offset with respect to the box
%                             centre. Data that are outliers are not
%                             included. The default is false.
%       showViolin          - If true, a 'violin' [1] kernel density
%                             outline of the data underlying each box will
%                             be plotted. The outline is calculated using
%                             the IOSR.STATISTICS.KERNELDENSITY function.
%                             The default is false.
%       style               - Determine whether to show additional x-axis 
%                             labels for the data. If set to 'hierarchy',
%                             additional hierarchical x-axis labels will be
%                             added for the plot. The labels can be set
%                             using the 'groupLabels' option. If set to
%                             'normal' (default) no additional labels will
%                             be plotted.
%       symbolColor         - Outlier marker color. The default is
%                             'auto' (see boxColor).
%       symbolMarker        - Marker used to denote outliers. The default
%                             is 'o'.
%       theme               - Specify a display theme to change multiple
%                             display properties. The options are:
%                                 'colorall'    : colored boxes (with 0.4 
%                                                 alpha), lines, and
%                                                 scatter markers (median
%                                                 line and mean marker are
%                                                 black)
%                                 'colorlines'  : clear boxes, black mean
%                                                 markers, colored lines
%                                                 and scatter markers
%                                 'colorboxes'  : colored boxes, black
%                                                 lines and markers, gray
%                                                 scatter markers
%                                 'default'     : clear boxes, colored
%                                                 median lines and markers,
%                                                 gray scatter markers
%                             Use the 'themeColors' property to specify the
%                             color(s) that are used.
%       themeColors           Colors used when creating the theme. The
%                             default is 'auto' (see boxColor).
%       violinBins          - If 'showViolin' is true, this specifes the
%                             bins used to calculate the kernel density.
%                             See IOSR.STATISTICS.KERNELDENSITY.
%       violinBinWidth      - If 'showViolin' is true, this specifes the
%                             bin width used to calculate the kernel
%                             density. See IOSR.STATISTICS.KERNELDENSITY.
%       violinColor         - Fill color for the violins. The default is
%                             'none' (see boxColor).
%       violinKernel        - If 'showViolin' is true, this specifes the
%                             kernel used to calculate the kernel density.
%                             See IOSR.STATISTICS.KERNELDENSITY.
%       xSeparator          - Logical value that when true adds a separator
%                             line between x groups. The default is false.
%       xSpacing            - Determine the x-axis spacing of boxes. By
%                             default ('x'), the data in x are used to
%                             determine the position of boxes on the x-axis
%                             (when x is numeric). Alternativley, when set
%                             to 'equal', boxes are equally-spaced, but the
%                             data in x are used to label the axis; the
%                             x-axis ticks are at 1:LENGTH(X).
%       handles             - Structure containing handles to the various
%                             objects that constitute the plot. The
%                             fields/handles are:
%                                 'axes'            : the parent axes of 
%                                                     the box plot
%                                 'fig'             : the parent figure of
%                                                     the box plot
%                                 'addPrctiles'     : chart line objects
%                                                     for each additional
%                                                     percentile marker
%                                 'addPrctilesTxt'  : text objects for each
%                                                     additional percentile
%                                                     marker
%                                 'box'             : patch objects for
%                                                     each box
%                                 'medianLines'     : line objects for each
%                                                     median line
%                                 'means'           : chart line objects
%                                                     for each mean marker
%                                 'notchLowerLines' : line objects for each
%                                                     lower notch line
%                                 'notchUpperLines' : line objects for each
%                                                     upper notch line
%                                 'upperWhiskers'   : line objects for each
%                                                     upper whisker line
%                                 'lowerWhiskers'   : line objects for each
%                                                     lower whisker line
%                                 'upperWhiskerTips': line objects for each
%                                                     upper whisker tip
%                                 'lowerWhiskerTips': line objects for each
%                                                     lower whisker line
%                                 'outliers'        : scatter objects for
%                                                     each set of outliers
%                                 'scatters'        : scatter objects for 
%                                                     each scatter overlay
%                                 'samplesTxt'      : text objects for each
%                                                     sample size display
%                                 'groupsTxt'       : text objects for each
%                                                     group label
%                                 'xseps'           : line objects for each
%                                                     x separator
%                                 'legend'          : the legend object
% 
%   These properties can be referenced using dot notation - e.g. H.BOXCOLOR
%   where H is an instance of the BOXPLOT object - or using the SET and GET
%   methods - e.g. GET(H,'BOXCOLOR'). Both methods are case-insensitive.
% 
%   Note that some handles will be empty unless the associated option is
%   turned on.
% 
%   Read-only properties:
%       x                   - The x data.
%       y                   - The y data.
%       weights             - Array giving the weights for the data in y.
%                             The array must be the same size as y. This
%                             option may be specified in the constructor.
%       statistics          - Structure containing the statistics used
%                             for the box plot. With the exception of
%                             'outliers', 'outliers_IX', 'addPrctiles', and
%                             'percentile' noted below, each field contains
%                             a 1-by-P-by-G-by-I-by-J... numeric array of
%                             values (identical to that returned by
%                             IOSR.STATISTICS.QUANTILE). The fields are:
%                                 'addPrctiles'   : an A-by-P-by-G... array
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
%                                 'PL'            : the lower percentiles
%                                                   specified by the
%                                                   'percentile' property,
%                                                   and the lower limit of
%                                                   the boxes (default is
%                                                   25th percentile)
%                                 'PU'            : the upper percentiles
%                                                   specified by the
%                                                   'percentile' property,
%                                                   and the upper limit of
%                                                   the boxes (default is
%                                                   75th percentile)
%                                 'Q1'            : the 25th percentiles
%                                 'Q3'            : the 75th percentiles
%                                 'IQR'           : the inter-quartile
%                                                   ranges
%                                 'mean'          : the mean values
%                                 'min'           : the minimum values
%                                                   (excl. outliers)
%                                 'max'           : the maximum values
%                                                   (excl. outliers)
%                                 'notch_u'       : the upper notch values
%                                 'notch_l'       : the lower notch values
%                                 'outliers'      : a 1-by-P-by-G cell
%                                                   array of outlier values
%                                 'outliers_IX'   : a logical array, the
%                                                   same size as Y, with
%                                                   true values indicating
%                                                   outliers
%                                 'std'          : the standard deviations
% 
%   In addition to the above specifications, some options can be specified
%   for each group. Parameters should be specified as a cell array of size
%   G-by-I-by-J... . These options are: 'boxColor', 'lineColor',
%   'lineStyle', 'meanColor', 'meanMarker, 'medianColor', 'notchLineColor',
%   'notchLineStyle', 'scatterMarker', 'scatterColor', 'symbolColor',
%   'symbolMarker', 'themeColor', and 'violinColor'.
% 
%   As noted above, colors may be specified as a colormap function handle
%   (e.g. @gray for grayscale colors). The function handle can refer to one
%   of the built-in colormap functions, or any other function capable of
%   generating valid colormaps. If the specified colormap is one of the
%   built-in colormaps 'hot' (@hot), 'gray' (@gray), 'bone' (@bone),
%   'copper' (@copper), or 'pink' (@pink), then boxPlot will take steps to
%   try to restrict the range of colors. Specifically, for fills (e.g.
%   'boxColor'), the colors will be restricted such that the minimum
%   luminance is 0.33 (in order to not mask dark lines). For lines, the
%   colors will be restricted such that the maximum luminance is 0.66 (in
%   order that lines are not masked by the white background).
%   
%   IOSR.STATISTICS.BOXPLOT methods:
%       boxPlot         - Create the box plot.
% 
%   See also IOSR.STATISTICS.TAB2BOX, IOSR.STATISTICS.QUANTILE, COLORMAP,
%       IOSR.STATISTICS.FUNCTIONALSPREADPLOT,
%       IOSR.STATISTICS.FUNCTIONALBOXPLOT, IOSR.STATISTICS.KERNELDENSITY.
% 
%   References
%   
%   [1] Hintze, Jerry L.; Nelson, Ray D. (1998). "Violin Plots: A Box
%       Plot-Density Trace Synergism". The American Statistician. 52 (2):
%       181?4.

%   Copyright 2016 University of Surrey.
    
    properties (AbortSet)
       addPrctiles = []             % Additional percentiles to plot.
       addPrctilesColors            % Colors for the additional percentile markers.
       addPrctilesLabels = {}       % Labels for additional percentiles.
       addPrctilesMarkers = {}      % Markers for additional percentiles.
       addPrctilesSize              % Size of the additional percentile markers.
       addPrctilesTxtSize           % Size of the additional percentile labels text.
       boxAlpha = 1                 % The transparency of the boxes.
       boxColor = 'none'            % Fill color of the boxes.
       boxWidth = 'auto'            % The width of the boxes.
       groupLabelFontSize = 9       % The font size of the group labels.
       groupLabelHeight = 'auto'    % The height of the area reserved for group labels.
       groupLabels = []             % Labels used to label the boxes for each x-group.
       groupWidth = 0.75            % The proportion of the x-axis interval across which each x-group of boxes should be spread.
       limit = '1.5IQR'             % Mode indicating the limits that define outliers.
       lineColor = 'k'              % Color of the box outlines and whiskers.
       lineStyle = '-'              % Style of the whisker lines.
       lineWidth = 1                % Width, in points, of the box outline, whisker lines, notch line, and outlier marker edges.
       meanColor = 'auto'           % Color of the mean marker.
       meanMarker = '+'             % Marker used for the mean.
       meanSize = 6                 % Size, in point units, of the mean markers.
       medianColor = 'auto'         % Color of the median line.
       method = 'R-8'               % The method used to calculate the quantiles.
       notch = false                % Whether the box should have a notch.
       notchDepth = 0.4             % Depth of the notch as a proportion of half the box width.
       notchLineColor = 'k'         % Color of the notch line.
       notchLineStyle = ':'         % Line style of the notch line.
       notchLine = false            % Whether to draw a horizontal line in the box at the extremes of the notch.
       outlierSize = 6^2            % Size, in square points, of the outlier markers.
       percentile = [25,75]         % Percentile limits of the box.
       sampleFontSize = 9           % Specify the font size of the sample size display.
       sampleSize = false           % Whether to display the sample size in the top-right of each box.
       scaleWidth = false           % Scale the width of each box according to the square root of the sample size.
       scatterAlpha = 1             % The transparency of the scatter markers.
       scatterColor = [.5 .5 .5]    % Scatter marker color for scatter plots of underlying data.
       scatterLayer = 'top'         % The layer of scatter plots with respect to the boxes.
       scatterMarker = 'x'          % Marker used for scatter plots of underlying data.
       scatterSize = 6^2            % Size, in square points, of the scatter markers.
       showLegend = false           % Draw a legend.
       showMean = false             % Display the mean of the data for each box.
       showOutliers = true          % Display outliers.
       showScatter = false          % Display a scatter plot of the underlying data for each box.
       showViolin = false           % Display a violin (kernel density) plot for the underlying data.
       style = 'normal'             % Determine whether to show additional x-axis labels for the data.
       symbolColor = 'auto'         % Outlier marker color.
       symbolMarker = 'o'           % Marker used to denote outliers.
       theme = 'default'            % Control a range of display properties.
       themeColors = 'auto'         % Colors used when creating the theme.
       violinBins = 'auto'          % The bins used to calculate the violins.
       violinBinWidth = 'auto'      % The width of the bins used to calculate the violins.
       violinAlpha = 1              % Alpha of the violins.
       violinColor = 'none'         % Color of the violins.
       violinWidth = 'auto'         % The width of the violins.
       violinKernel = 'normal'      % The violin kernel used to calculate the kernel density.
       xSeparator = false           % Add a separator line between x groups.
       xSpacing = 'x'               % Determine the x-axis spacing of boxes.
    end
    
    properties (Access = private)
        diffx                       % smallest difference of x axis data
        groupXticks                 % x ticks for every box
        groupDims                   % dimensions of data for each x group
        groupRange                  % scale value range of each x group
        xlabels                     % labels for the x axis
        xticks                      % x-axis tick points
    end
    
    properties (SetAccess = private, Dependent, Hidden)
        addPrctilesHandle = []      % Chart line objects objects for the additional percentiles markers.
        addPrctilesTxtHandle = []   % Text objects for each additional percentile label.
        boxHandles = []             % Patch objects for each box.
        mLineHandles = []           % Line objects for each median line.
        meanHandles = []            % Chart line objects for each mean marker.
        notchLowerLineHandles = []  % Line objects for each lower notch line.
        notchUpperLineHandles = []  % Line objects for each upper notch line.
        upperWhiskerHandles = []    % Line objects for each upper whisker line.
        lowerWhiskerHandles = []    % Line objects for each lower whisker line.
        upperWhiskerTipHandles = [] % Line objects for each upper whisker tip.
        lowerWhiskerTipHandles = [] % Line objects for each lower whisker tip.
        outliersHandles = []        % Scatter objects for each set of outliers.
        scatterHandles = []         % Scatter objects for each scatter overlay.
        sampleTxtHandes = []        % Text objects for each sample size display.
        groupTxtHandes = []         % Text objects for each group label.
        xsepHandles = []            % Line objects for each x separator.
        legendHandle = []           % Legend object.
    end
    
    methods
        
        function obj = boxPlot(varargin)
        % BOXPLOT Draw a box plot.
        %
        %   IOSR.STATISTICS.BOXPLOT(Y) produces a box plot of the data in
        %   Y. If Y is a vector, there is just one box. If Y is a matrix,
        %   there is one box per column. If Y is an
        %   N-by-P-by-G-by-I-by-J... array then G*I*J*... boxes are plotted
        %   hierarchically for each column P; i.e. J boxes are plotted for
        %   each index of I, and J*I boxes are plotted for each index of G.
        %   For each box, the central mark is the median of the
        %   data/column, the edges of the box are the 25th and 75th
        %   percentiles, the whiskers extend to the most extreme data
        %   points not considered outliers, and outliers are plotted
        %   individually. NaNs are excluded from the data.
        %
        %   Tabular data can be arranged into the appropriate format using
        %   the IOSR.STATISTICS.TAB2BOX function.
        %
        %   IOSR.STATISTICS.BOXPLOT(X,Y) specifies the x-axis values for
        %   each box. X should be a vector, with as many elements as Y has
        %   columns. The default is 1:SIZE(Y,2).
        %
        %   IOSR.STATISTICS.BOXPLOT(...,'PARAMETER',VALUE) allows the
        %   plotting options to be specified when the plot is constructed.
        %
        %   IOSR.STATISTICS.BOXPLOT(AX,...) creates the box plot in the
        %   axes specified by AX.
        %
        %   Examples
        % 
        %     Example 1: Basic grouped box plot with legend
        % 
        %       y = randn(50,3,3);
        %       x = [1 2 3.5];
        %       y(1:25) = NaN;
        % 
        %       figure;
        %       h = iosr.statistics.boxPlot(x,y,...
        %           'symbolColor','k',...
        %           'medianColor','k',...
        %           'symbolMarker',{'+','o','d'},...
        %           'boxcolor',{[1 0 0]; [0 1 0]; [0 0 1]},...
        %           'groupLabels',{'y1','y2','y3'},...
        %           'showLegend',true);
        %       box on
        % 
        %     Example 2: Grouped box plot with overlayed data
        % 
        %       figure;
        %       iosr.statistics.boxPlot(x,y,...
        %           'symbolColor','k',...
        %           'medianColor','k',...
        %           'symbolMarker',{'+','o','d'},...
        %           'boxcolor','auto',...
        %           'showScatter',true);
        %       box on
        % 
        %     Example 3: Grouped box plot with displayed sample sizes
        %       and variable widths
        % 
        %       figure;
        %       iosr.statistics.boxPlot(x,y,...
        %           'medianColor','k',...
        %           'symbolMarker',{'+','o','d'},...
        %           'boxcolor','auto',...
        %           'sampleSize',true,...
        %           'scaleWidth',true);
        %       box on
        % 
        %     Example 4: Grouped notched box plot with x separators and
        %       hierarchical labels
        % 
        %       figure;
        %       iosr.statistics.boxPlot({'A','B','C'},y,...
        %           'notch',true,...
        %           'medianColor','k',...
        %           'symbolMarker',{'+','o','d'},...
        %           'boxcolor','auto',...
        %           'style','hierarchy',...
        %           'xSeparator',true,...
        %           'groupLabels',{{'Group 1','Group 2','Group 3'}});
        %       box on
        %
        %     Example 5: Box plot with legend labels from data
        % 
        %       % load data
        %       % (requires Statistics or Machine Learning Toolbox)
        %       load carbig
        % 
        %       % arrange data
        %       [y,x,g] = iosr.statistics.tab2box(Cylinders,MPG,when);
        %   
        %       % sort
        %       IX = [1 3 2]; % order
        %       g = g{1}(IX);
        %       y = y(:,:,IX);
        % 
        %       % plot
        %       figure
        %       h = iosr.statistics.boxPlot(x,y,...
        %           'symbolColor','k','medianColor','k','symbolMarker','+',...
        %           'boxcolor',{[1 1 1],[.75 .75 .75],[.5 .5 .5]},...
        %           'scalewidth',true,'xseparator',true,...
        %           'groupLabels',g,'showLegend',true);
        %       box on
        %       title('MPG by number of cylinders and period')
        %       xlabel('Number of cylinders')
        %       ylabel('MPG')
        %
        %     Example 6: Box plot calculated from weighted quantiles
        %
        %       % load data
        %       load carbig
        %       
        %       % random weights
        %       weights = rand(size(MPG));
        %       
        %       % arrange data
        %       [y,x,g] = iosr.statistics.tab2box(Cylinders,MPG,when);
        %       weights_boxed = iosr.statistics.tab2box(Cylinders,weights,when);
        %       
        %       % plot
        %       figure
        %       h = iosr.statistics.boxPlot(x,y,'weights',weights_boxed);
        %
        %     Example 7: Draw a violin plot
        %       y = randn(50,3,3);
        %       x = [1 2 3.5];
        %       y(1:25) = NaN;
        %       figure('color','w');
        %       h2 = iosr.statistics.boxPlot(x,y, 'showViolin', true, 'boxWidth', 0.025, 'showOutliers', false);
        %       box on
            
            if nargin > 0
        
                %% check for deprecated properties

                % check for scatter option
                scatterIX = strcmpi('scatter',varargin);
                if any(scatterIX)
                    warning('The ''scatter'' property is deprecated. Use ''showScatter'' instead.');
                    varargin{scatterIX} = 'showScatter';
                end

                % check for mean option
                meanIX = strcmpi('mean',varargin);
                if any(meanIX)
                    warning('The ''mean'' property is deprecated. Use ''showMean'' instead.');
                    varargin{meanIX} = 'showMean';
                end

                %% set x, y, and dims

                % check for input data
                start = obj.getXY(varargin{:});
                obj.groupDims = obj.ydims(3:end);

                % create initial theme
                obj.createTheme();

                % check weights
                obj.checkWeights();

                % set properties from varargin
                obj.setProperties(start,nargin,varargin);

                % remove NaN columns
                obj.removeNaN();

                % set x axis properties and labels
                obj.setXprops();

                %% statistics

                % calculate statistics
                obj.calculateStats();

                %% draw

                % set handles
                obj.parseAxesHandle(varargin{:});

                % draw the box plot
                obj.draw('all');

                % add listener to ylim so that some properties can be redrawn
                addlistener(handle(obj.handles.axes),'YLim','PostSet',...
                    @(hProp,eventData) obj.eventHandler(hProp,eventData));
                
            end
            
        end
        
        %% dependent handle getters
        
        % these are all deprecated and the handles moved to the obj.handles struct
        
        function val = get.addPrctilesHandle(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.addPrctilesHandle'' is deprecated. Use ''obj.handles.addPrctiles'' instead.');
            val = obj.handles.addPrctiles;
        end
        
        function val = get.addPrctilesTxtHandle(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.addPrctilesTxtHandle'' is deprecated. Use ''obj.handles.addPrctilesTxt'' instead.');
            val = obj.handles.addPrctilesTxt;
        end
        
        function val = get.boxHandles(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.boxHandles'' is deprecated. Use ''obj.handles.box'' instead.');
            val = obj.handles.box;
        end
        
        function val = get.mLineHandles(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.mLineHandles'' is deprecated. Use ''obj.handles.medianLines'' instead.');
            val = obj.handles.medianLines;
        end
        
        function val = get.meanHandles(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.meanHandles'' is deprecated. Use ''obj.handles.means'' instead.');
            val = obj.handles.means;
        end
        
        function val = get.notchLowerLineHandles(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.notchLowerLineHandles'' is deprecated. Use ''obj.handles.notchLowerLines'' instead.');
            val = obj.handles.notchLowerLines;
        end
        
        function val = get.notchUpperLineHandles(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.notchUpperLineHandles'' is deprecated. Use ''obj.handles.notchUpperLines'' instead.');
            val = obj.handles.notchUpperLines;
        end
        
        function val = get.upperWhiskerHandles(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.upperWhiskerHandles'' is deprecated. Use ''obj.handles.upperWhiskers'' instead.');
            val = obj.handles.upperWhiskers;
        end
        
        function val = get.lowerWhiskerHandles(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.lowerWhiskerHandles'' is deprecated. Use ''obj.handles.lowerWhiskers'' instead.');
            val = obj.handles.lowerWhiskers;
        end
        
        function val = get.upperWhiskerTipHandles(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.upperWhiskerTipHandles'' is deprecated. Use ''obj.handles.upperWhiskerTips'' instead.');
            val = obj.handles.upperWhiskerTips;
        end
        
        function val = get.lowerWhiskerTipHandles(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.lowerWhiskerTipHandles'' is deprecated. Use ''obj.handles.lowerWhiskerTips'' instead.');
            val = obj.handles.lowerWhiskerTips;
        end
        
        function val = get.outliersHandles(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.outliersHandles'' is deprecated. Use ''obj.handles.outliers'' instead.');
            val = obj.handles.outliers;
        end
        
        function val = get.scatterHandles(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.scatterHandles'' is deprecated. Use ''obj.handles.scatters'' instead.');
            val = obj.handles.scatters;
        end
        
        function val = get.sampleTxtHandes(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.sampleTxtHandes'' is deprecated. Use ''obj.handles.samplesTxt'' instead.');
            val = obj.handles.samplesTxt;
        end
        
        function val = get.groupTxtHandes(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.groupTxtHandes'' is deprecated. Use ''obj.handles.groupsTxt'' instead.');
            val = obj.handles.groupsTxt;
        end
        
        function val = get.xsepHandles(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.xsepHandles'' is deprecated. Use ''obj.handles.xseps'' instead.');
            val = obj.handles.xseps;
        end
        
        function val = get.legendHandle(obj)
            warning('iosr:boxPlot:deprecatedProp','''obj.legendHandle'' is deprecated. Use ''obj.handles.legend'' instead.');
            val = obj.handles.legend;
        end
        
        %% accessor functions
        
        % set/get additional percentiles
        
        function val = get.addPrctiles(obj)
            val = obj.addPrctiles;
        end
        
        function set.addPrctiles(obj,val)
            if ~isempty(val)
                assert(isvector(val) && isnumeric(val), 'iosr:boxPlot:addPrctilesVector', '''ADDPRCTILES'' should be a numeric vector')
                assert(all(val<=100) && all(val>=0), 'iosr:boxPlot:addPrctilesRange', '''ADDPRCTILES'' values should be in the interval [0,100].')
            end
            obj.addPrctiles = val;
            obj.calculateStats();
            obj.draw('addPrctiles');
        end
        
        % set/get additional percentiles colors
        
        function val = get.addPrctilesColors(obj)
            val = obj.checkAddPrctileProp(obj.addPrctilesColors,{'k'});
        end
        
        function set.addPrctilesColors(obj,val)
            if ~isempty(val)
                assert(iscell(val), 'iosr:boxPlot:addPrctilesColors', '''ADDPRCTILESCOLORS'' should be a cell array')
            end
            obj.addPrctilesColors = val;
            obj.draw();
        end
        
        % set/get additional percentiles labels
        
        function val = get.addPrctilesLabels(obj)
            val = obj.checkAddPrctileProp(obj.addPrctilesLabels,{''});
        end
        
        function set.addPrctilesLabels(obj,val)
            if ~isempty(val)
                assert(iscellstr(val), 'iosr:boxPlot:addPrctilesLabels', '''ADDPRCTILESLABELS'' should be a cell array of strings')
            end
            obj.addPrctilesLabels = val;
            obj.draw();
        end
        
        % set/get additional percentiles markers
        
        function val = get.addPrctilesMarkers(obj)
            val = obj.checkAddPrctileProp(obj.addPrctilesMarkers,{'*'});
        end
        
        function set.addPrctilesMarkers(obj,val)
            if ~isempty(val)
                assert(iscellstr(val), 'iosr:boxPlot:addPrctilesMarkers', '''ADDPRCTILESMARKERS'' should be a cell array of strings')
            end
            obj.addPrctilesMarkers = val;
            obj.draw();
        end
        
        % set/get additional percentiles marker size
        
        function val = get.addPrctilesSize(obj)
            val = obj.checkAddPrctileProp(obj.addPrctilesSize,6);
        end
        
        function set.addPrctilesSize(obj,val)
            if ~isempty(val)
                assert(isnumeric(val) && isvector(val), 'iosr:boxPlot:addPrctilesSize', '''ADDPRCTILESSIZE'' should be a numeric vector')
            end
            obj.addPrctilesSize = val;
            obj.draw();
        end
        
        % set additional percentiles label font size
        
        function set.addPrctilesTxtSize(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:boxPlot:addPrctilesTxtSize', '''ADDPRCTILESTXTSIZE'' should be a numeric scalar')
            obj.addPrctilesTxtSize = val;
            obj.draw();
        end
        
        % set box alpha
        
        function set.boxAlpha(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:boxPlot:boxAlpha', '''BOXALPHA'' must be a numeric scalar')
            obj.boxAlpha = val;
            obj.draw('legend');
        end
        
        % set/get box color
        
        function val = get.boxColor(obj)
            val = obj.checkColor(obj.boxColor,'fill');
            val = obj.groupOption(val,'boxColor');
        end
        
        function set.boxColor(obj,val)
            obj.boxColor = val;
            obj.draw('legend');
        end
        
        % set box width
        
        function set.boxWidth(obj,val)
            assert((isnumeric(val) && isscalar(val)) || strcmpi(val,'auto'), 'iosr:boxPlot:boxWidth', '''BOXWIDTH'' must be a numeric scalar or ''auto''.')
            obj.boxWidth = val;
            obj.draw('all');
        end
        
        % set/get group labels
        
        function val = get.groupLabels(obj)
            if prod(obj.groupDims)==numel(obj.groupLabels) && numel(obj.groupLabels)>1
                val = {obj.groupLabels};
            else
                val = obj.groupLabels;
            end
            if ~isempty(val) % use input
                assert(isvector(val) && iscell(val),...
                    'iosr:boxPlot:groupLabelsType', ...
                    'The GROUPLABELS option should be a cell vector');
                valSize = cellfun(@length,val);
                assert(prod(obj.outDims(3:end))==prod(valSize), ...
                    'iosr:boxPlot:groupLabelsSize', ...
                    ['The GROUPLABELS option should be a cell vector; ' ...
                    'the Nth element should contain a vector of length SIZE(Y,N+2)'])
                assert(isequal(obj.outDims(3:end),valSize(1:length(obj.outDims)-2)), ...
                    'iosr:boxPlot:groupLabelsSize', ...
                    ['The GROUPLABELS option should be a cell vector; ' ...
                    'the Nth element should contain a vector of length SIZE(Y,N+2)'])
            else % create placeholder labels
                val = cell(1,length(obj.groupDims));
                c = 1;
                for m = 1:length(val)
                    val{m} = cell(1,obj.groupDims(m));
                    for n = 1:obj.groupDims(m)
                        val{m}{n} = sprintf('Data %d',c);
                        c = c+1;
                    end
                end
            end
        end
        
        function set.groupLabels(obj,val)
            obj.groupLabels = val;
            obj.draw('style','legend');
        end
        
        % set group label font size
        
        function set.groupLabelFontSize(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:boxPlot:groupLabelsFontSize', '''GROUPLABELFONTSIZE'' must be a numeric scalar.')
            obj.groupLabelFontSize = val;
            obj.draw();
        end
        
        % set group label height
        
        function set.groupLabelHeight(obj,val)
            assert(strcmp(val,'auto') || (isnumeric(val) && isscalar(val)), 'iosr:boxPlot:groupLabelHeight', '''GROUPLABELHEIGHT'' must be ''auto'' or a numeric scalar.')
            obj.groupLabelHeight = val;
            obj.draw('style');
        end
        
        % set group width
        
        function set.groupWidth(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:boxPlot:groupWidth', '''GROUPWIDTH'' must be a numeric scalar')
            obj.groupWidth = val;
            obj.setXprops();
            obj.draw('all');
        end
        
        % set stats limit
        
        function set.limit(obj,val)
            if isnumeric(val)
                assert(isnumeric(val) && numel(val)==2, 'iosr:boxPlot:limitSize', '''LIMIT'' must be a two-element numeric vector.')
                assert(all(val<=100) && all(val>=0), 'iosr:boxPlot:limitRange', '''LIMIT'' values should be in the interval [0,100].')
                if all(val<1)
                    warning('iosr:boxPlot:limitRange','''LIMIT'' values should be in the interval [0,100].')
                end
            else
                assert(ischar(val), 'iosr:boxPlot:limitType', '''LIMIT'' must be a char or numeric array.')
                assert(any(strcmpi(val,{'1.5IQR','3IQR','none'})), 'iosr:boxPlot:limitUnknownOption', '''LIMIT'' parameter not recognised.')
            end
            obj.limit = val;
            obj.calculateStats();
            obj.draw('whiskers','outliers','scatter');
        end
        
        % get/set line colors
        
        function val = get.lineColor(obj)
            val = obj.checkColor(obj.lineColor,'line');
            val = obj.groupOption(val,'lineColor');
        end
        
        function set.lineColor(obj,val)
            obj.lineColor = val;
            obj.draw('legend');
        end
        
        % get/set line style
        
        function val = get.lineStyle(obj)
            val = obj.groupOption(obj.lineStyle,'lineStyle');
        end
        
        function set.lineStyle(obj,val)
            assert(ischar(val) || iscellstr(val), 'iosr:boxPlot:lineStyle', '''LINESTYLE'' must be a char array or cell array of strings.')
            obj.lineStyle = val;
            obj.draw('legend');
        end
        
        % set line width
        
        function set.lineWidth(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:boxPlot:lineWidth', '''LINEWIDTH'' must be a numeric scalar.')
            obj.lineWidth = val;
            obj.draw('all');
        end
        
        % get/set mean color
        
        function val = get.meanColor(obj)
            val = obj.checkColor(obj.meanColor,'line');
            val = obj.groupOption(val,'meanColor');
        end
        
        function set.meanColor(obj,val)
            obj.meanColor = val;
            obj.draw();
        end
        
        % get/set mean marker
        
        function val = get.meanMarker(obj)
            val = obj.groupOption(obj.meanMarker,'meanMarker');
        end
        
        function set.meanMarker(obj,val)
            assert(ischar(val) || iscellstr(val), 'iosr:boxPlot:meanMarker', '''MEANMARKER'' must be a char array or cell array of strings.')
            obj.meanMarker = val;
            obj.draw();
        end
        
        % set mean size
        
        function set.meanSize(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:boxPlot:meanSize', '''MEANSIZE'' must be a numeric scalar.')
            obj.meanSize = val;
            obj.draw();
        end
        
        % get/set median color
        
        function val = get.medianColor(obj)
            val = obj.checkColor(obj.medianColor,'line');
            val = obj.groupOption(val,'medianColor');
        end
        
        function set.medianColor(obj,val)
            obj.medianColor = val;
            obj.draw('legend');
        end
        
        % set stats method
        
        function set.method(obj,val)
            assert(ischar(val), 'iosr:boxPlot:method', '''METHOD'' must be a char array.')
            obj.method = val;
            obj.calculateStats();
            obj.draw('all');
        end
        
        % set notch
        
        function set.notch(obj,val)
            assert(islogical(val) && numel(val)==1, 'iosr:boxPlot:notch', '''NOTCH'' must be logical.')
            obj.notch = val;
            obj.draw('boxes');
        end
        
        % set notch depth
        
        function set.notchDepth(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:boxPlot:notchDepth', '''NOTCHDEPTH'' must be a numeric scalar')
            obj.notchDepth = val;
            obj.draw('boxes');
        end
        
        % get/set notch line color
        
        function val = get.notchLineColor(obj)
            val = obj.checkColor(obj.notchLineColor,'line');
            val = obj.groupOption(val,'notchLineColor');
        end
        
        function set.notchLineColor(obj,val)
            obj.notchLineColor = val;
            obj.draw();
        end
        
        % get/set notch line style
        
        function val = get.notchLineStyle(obj)
            val = obj.groupOption(obj.notchLineStyle,'notchLineStyle');
        end
        
        function set.notchLineStyle(obj,val)
            assert(ischar(val) || iscellstr(val), 'iosr:boxPlot:notchLineStyle', '''NOTCHLINESTYLE'' must be a char array or cell array of strings.')
            obj.notchLineStyle = val;
            obj.draw();
        end
        
        % set notch line
        
        function set.notchLine(obj,val)
            assert(islogical(val) && numel(val)==1, 'iosr:boxPlot:notchLine', '''NOTCHLINE'' must be logical.')
            obj.notchLine = val;
            obj.draw('boxes');
        end
        
        % set outlier size
        
        function set.outlierSize(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:boxPlot:outlierSize', '''OUTLIERSIZE'' must be a numeric scalar.')
            obj.outlierSize = val;
            obj.draw();
        end
        
        % set percentile
        
        function set.percentile(obj,val)
            assert(isnumeric(val) && numel(val)==2, 'iosr:boxPlot:percentileSize', '''PERCENTILE'' must be a two-element numeric vector.')
            assert(all(val<=100) && all(val>=0), 'iosr:boxPlot:percentileRange', '''PERCENTILE'' values should be in the interval [0,100].')
            if all(val<1)
                warning('iosr:boxPlot:percentileRange', '''PERCENTILE'' values should be in the interval [0,100].')
            end
            obj.percentile = val;
            obj.calculateStats();
            obj.draw('all');
        end
        
        % set sample font size
        
        function set.sampleFontSize(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:boxPlot:sampleFontSize', '''SAMPLEFONTSIZE'' must be a numeric scalar.')
            obj.sampleFontSize = val;
            obj.draw();
        end
        
        % set sample size option
        
        function set.sampleSize(obj,val)
            assert(islogical(val) && numel(val)==1, 'iosr:boxPlot:sampleSize', '''SAMPLESIZE'' must be logical.')
            obj.sampleSize = val;
            obj.draw('boxes');
        end
        
        % set scale width option
        
        function set.scaleWidth(obj,val)
            assert(islogical(val) && numel(val)==1, 'iosr:boxPlot:scaleWidth', '''SCALEWIDTH'' must be logical.')
            obj.scaleWidth = val;
            obj.draw('boxes','scatter','outliers');
        end
        
        % set legend
        
        function set.showLegend(obj,val)
            assert(islogical(val) && numel(val)==1, 'iosr:boxPlot:showLegend', '''SHOWLEGEND'' must be logical.')
            obj.showLegend = val;
            obj.draw('legend');
        end
        
        % set show mean option
        
        function set.showMean(obj,val)
            assert(islogical(val) && numel(val)==1, 'iosr:boxPlot:showMean', '''SHOWMEAN'' must be logical.')
            obj.showMean = val;
            obj.draw('means');
        end
        
        % set show outliers option
        
        function set.showOutliers(obj,val)
            assert(islogical(val) && numel(val)==1, 'iosr:boxPlot:showOutliers', '''SHOWOUTLIERS'' must be logical.')
            obj.showOutliers = val;
            obj.draw('outliers');
        end
        
        % set show scatter option
        
        function set.showScatter(obj,val)
            assert(islogical(val) && numel(val)==1, 'iosr:boxPlot:showScatter', '''SHOWSCATTER'' must be logical.')
            obj.showScatter = val;
            obj.draw('scatter');
        end
        
        % set scatter alpha
        
        function set.scatterAlpha(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:boxPlot:scatterAlpha', '''SCATTERALPHA'' must be a numeric scalar')
            obj.scatterAlpha = val;
            obj.draw();
        end
        
        % get/set scatter color
        
        function val = get.scatterColor(obj)
            val = obj.checkColor(obj.scatterColor,'line');
            val = obj.groupOption(val,'scatterColor');
        end
        
        function set.scatterColor(obj,val)
            obj.scatterColor = val;
            obj.draw();
        end
        
        % set scatter layer
        
        function set.scatterLayer(obj,val)
            assert(ischar(val), 'iosr:boxPlot:scatterLayer', '''SCATTERLAYER'' must be a char array.')
            obj.scatterLayer = val;
            obj.draw();
        end
        
        % get/set scatter marker
        
        function val = get.scatterMarker(obj)
            val = obj.groupOption(obj.scatterMarker,'scatterMarker');
        end
        
        function set.scatterMarker(obj,val)
            assert(ischar(val) || iscellstr(val), 'iosr:boxPlot:scatterMarker', '''SCATTERMARKER'' must be a char array or cell array of strings.')
            obj.scatterMarker = val;
            obj.draw();
        end
        
        % set scatter marker size
        
        function set.scatterSize(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:boxPlot:scatterSize', '''SCATTERSIZE'' must be a numeric scalar.')
            obj.scatterSize = val;
            obj.draw();
        end
        
        % violin
        
        function set.showViolin(obj,val)
            assert(islogical(val) && isscalar(val), 'iosr:boxPlot:showViolin', '''VIOLIN'' must be logical.')
            obj.showViolin = val;
            obj.draw('violin','whiskers');
        end
        
        % set style
        
        function set.style(obj,val)
            assert(ischar(val), 'iosr:boxPlot:styleType', '''STYLE'' must be a char array.')
            assert(any(strcmpi(val,{'normal','hierarchy'})), 'iosr:boxPlot:styleOption', '''STYLE'' must be either ''normal'' or ''hierarchy''.')
            obj.style = val;
            obj.draw('style');
        end
        
        % get/set outlier symbol color
        
        function val = get.symbolColor(obj)
            val = obj.checkColor(obj.symbolColor,'line');
            val = obj.groupOption(val,'symbolColor');
        end
        
        function set.symbolColor(obj,val)
            obj.symbolColor = val;
            obj.draw();
        end
        
        % get/set outlier symbol marker
        
        function val = get.symbolMarker(obj)
            val = obj.groupOption(obj.symbolMarker,'symbolMarker');
        end
        
        function set.symbolMarker(obj,val)
            assert(ischar(val) || iscellstr(val), 'iosr:boxPlot:symbolMarker', '''SYMBOLMARKER'' must be a char array or cell array of strings.')
            obj.symbolMarker = val;
            obj.draw();
        end
        
        % set theme
        
        function set.theme(obj,val)
            assert(ischar(val), 'iosr:boxPlot:theme', '''THEME'' must be a char array.')
            obj.theme = val;
            obj.createTheme();
            obj.draw('legend');
        end
        
        % theme color
        
        function val = get.themeColors(obj)
            switch lower(obj.theme)
                case {'graylines', 'colorlines', 'default'}
                    val = obj.checkColor(obj.themeColors,'lines');
                case {'grayboxes', 'colorall', 'colorboxes'}
                    val = obj.checkColor(obj.themeColors,'fill');
            end
            val = obj.groupOption(val,'themeColors');
        end
        
        function set.themeColors(obj,val)
            obj.themeColors = val;
            obj.createTheme();
            obj.draw('all');
        end
        
        % set violin props
        
        function set.violinBins(obj,val)
            assert(isnumeric(val) || strcmpi(val,'auto'), 'iosr:boxPlot:violinBins', '''VIOLINBINS'' must be numeric or ''auto''.')
            obj.violinBins = val;
            obj.draw('violin','whiskers');
        end
        
        function set.violinBinWidth(obj,val)
            assert((isnumeric(val) && isscalar(val)) || strcmpi(val,'auto'), 'iosr:boxPlot:violinBinWidth', '''VIOLINBINWIDTH'' must be a numeric scalar or ''auto''.')
            obj.violinBinWidth = val;
            obj.draw('violin','whiskers');
        end
        
        function set.violinKernel(obj,val)
            assert(ischar(val), 'iosr:boxPlot:violinKernel', '''VIOLINKERNEL'' must be a string.')
            obj.violinKernel = val;
            obj.draw('violin','whiskers');
        end
        
        function set.violinAlpha(obj,val)
            assert(isnumeric(val) && isscalar(val), 'iosr:boxPlot:violinAlpha', '''VIOLINALPHA'' must be a numeric scalar')
            obj.violinAlpha = val;
            obj.draw('legend');
        end
        
        function val = get.violinColor(obj)
            val = obj.checkColor(obj.violinColor,'fill');
            val = obj.groupOption(val,'violinColor');
        end
        
        function set.violinColor(obj,val)
            obj.violinColor = val;
            obj.draw('legend');
        end
        
        function set.violinWidth(obj,val)
            assert((isnumeric(val) && isscalar(val)) || strcmpi(val,'auto'), 'iosr:boxPlot:violinWidth', '''VIOLINWIDTH'' must be a numeric scalar or ''auto''.')
            obj.violinWidth = val;
            obj.draw('violin','whiskers');
        end
        
        % set x separator option
        
        function set.xSeparator(obj,val)
            assert(islogical(val) && isscalar(val), 'iosr:boxPlot:xSeparator', '''XSEPARATOR'' must be logical.')
            obj.xSeparator = val;
            obj.draw('xsep');
        end
        
        % set x spacing option
        
        function set.xSpacing(obj,val)
            assert(any(strcmpi(val,{'x','equal'})), 'iosr:boxPlot:xSpacing', '''XSPACING'' parameter not recognised.')
            obj.xSpacing = val;
            obj.setXprops();
            obj.draw('all');
        end
        
        
    end
    
    methods (Access = protected)
        
        function draw(obj,varargin)
        %DRAW main draw function
            
            if isfield(obj.handles,'axes')
                subidx = cell(1,length(obj.outDims));
                hold on;
                for n = 1:prod(obj.outDims)
                    % get indices
                    [subidx{:}] = ind2sub(obj.outDims, n);
                    subidxAll = [{':'} subidx(2:end)];
                    gidx = subidx(3:end);

                    % determine where to put tick marks for group labels
                    obj.groupXticks(subidx{:}) = obj.xticks(subidx{2})+obj.getOffset(subidx);
                    
                    % plot individual components
                    if any(strcmpi('violin',varargin)) || any(strcmpi('all',varargin))
                        obj.drawViolin(subidx,subidxAll);
                    end
                    if any(strcmpi('boxes',varargin)) || any(strcmpi('all',varargin))
                        obj.drawBox(subidx,subidxAll);
                    end
                    if any(strcmpi('whiskers',varargin)) || any(strcmpi('all',varargin))
                        obj.drawWhiskers(subidx);
                    end
                    if any(strcmpi('outliers',varargin)) || any(strcmpi('all',varargin))
                        obj.drawOutliers(subidx);
                    end
                    if any(strcmpi('scatter',varargin)) || any(strcmpi('all',varargin))
                        obj.drawScatter(subidx,gidx,subidxAll);
                    end
                    if any(strcmpi('means',varargin)) || any(strcmpi('all',varargin))
                        obj.drawMean(subidx);
                    end
                    if any(strcmpi('addPrctiles',varargin)) || any(strcmpi('all',varargin))
                        obj.drawAddPrctiles(subidx);
                    end
                    % change graphics options
                    obj.drawGroupGraphics(subidx,gidx);
                end
                % style things
                if any(strcmpi('xsep',varargin)) || any(strcmpi('all',varargin))
                    obj.drawXsepatator();
                end
                if any(strcmpi('style',varargin)) || any(strcmpi('all',varargin))
                    obj.drawStyle();
                end
                if any(strcmpi('legend',varargin)) || any(strcmpi('all',varargin))
                    obj.drawLegend();
                end
                % change graphics options
                obj.drawGlobalGraphics();
                % set layering
                if isfield(obj.handles,'scatters')
                    uistack(obj.handles.scatters(:),obj.scatterLayer);
                end
                if isfield(obj.handles,'samplesTxt')
                    uistack(obj.handles.samplesTxt(:),'top');
                end
                if isfield(obj.handles,'means')
                    uistack(obj.handles.means(:),'top');
                end
                if isfield(obj.handles,'addPrctiles')
                    uistack(obj.handles.addPrctiles(:),'top');
                end
                if isfield(obj.handles,'addPrctilesTxt')
                    uistack(obj.handles.addPrctilesTxt(:),'top');
                end
                hold off;
            end
        end
        
        function drawViolin(obj,subidx,subidxAll)
        %DRAWVIOLIN Draw the violin
        
            try % to delete the handles first
                delete(obj.handles.violin(subidx{:}));
            catch
            end
            
            if obj.showViolin
                
                vBinWidth = obj.replaceAuto(...
                    obj.violinBinWidth, ...
                    [], ...
                    'Unknown ''VIOLINBINWIDTH'' parameter.');
                
                vBins = obj.replaceAuto(...
                    obj.violinBins, ...
                    [], ...
                    'Unknown ''VIOLINBINS'' parameter.');

                [d, xd] = iosr.statistics.kernelDensity(obj.y(subidxAll{:}), vBins, vBinWidth, obj.violinKernel);
                halfviolinwidth = obj.calcHalfViolinWidth(subidx);
                d = halfviolinwidth .* (d ./ max(d));
                obj.handles.violin(subidx{:}) = patch(...
                    [d; -flipud(d)] + obj.xticks(subidx{2}) + obj.getOffset(subidx), ...
                    [xd; flipud(xd)],'w','Parent',obj.handles.axes);
                
            end
            
        end
        
        function drawBox(obj,subidx,subidxAll)
        %DRAWBOX draw boxes
            
            try % to delete the handles first
                delete(obj.handles.box(subidx{:}));
            catch
            end
            try
                delete(obj.handles.medianLines(subidx{:}));
            catch
            end
            try
                delete(obj.handles.notchUpperLines(subidx{:}));
            catch
            end
            try
                delete(obj.handles.notchLowerLines(subidx{:}));
            catch
            end
            try
                delete(obj.handles.samplesTxt(subidx{:}));
            catch
            end
            
            % determine half box width
            [~, halfboxwidth] = obj.calcBoxWidths(subidx);

            % notch depth
            notchdepth = obj.notchDepth*halfboxwidth;

            % main nodes
            X = obj.xticks(subidx{2}) + obj.getOffset(subidx);
            q1 = obj.statistics.PL(subidx{:});
            q3 = obj.statistics.PU(subidx{:});
            md = obj.statistics.median(subidx{:});
            nu = obj.statistics.notch_u(subidx{:});
            nl = obj.statistics.notch_l(subidx{:});
            
            % set nodes
            if obj.notch % if notch requested
                Xnodes = [X-halfboxwidth X-halfboxwidth X-halfboxwidth+notchdepth X-halfboxwidth X-halfboxwidth ...
                    X+halfboxwidth X+halfboxwidth X+halfboxwidth-notchdepth X+halfboxwidth X+halfboxwidth];
                Ynodes = [q1 nl md nu q3 q3 nu md nl q1];
            else % if not
                Xnodes = [X-halfboxwidth X-halfboxwidth X+halfboxwidth X+halfboxwidth];
                Ynodes = [q1 q3 q3 q1];
                
            end
            
            % draw box and set props
            obj.handles.box(subidx{:}) = patch(Xnodes,Ynodes,'w','Parent',obj.handles.axes);
            
            if obj.notchLine % if notchLine requested
                obj.handles.notchLowerLines(subidx{:}) = line([X-halfboxwidth X+halfboxwidth],[nl nl],...
                    'Parent',obj.handles.axes);
                obj.handles.notchUpperLines(subidx{:}) = line([X-halfboxwidth X+halfboxwidth],[nu nu],...
                    'Parent',obj.handles.axes);
            end
            
            if obj.notch % if notch requested
                % median
                obj.handles.medianLines(subidx{:}) = line([X-halfboxwidth+notchdepth X+halfboxwidth-notchdepth], [md md]);
            else
                % median
                obj.handles.medianLines(subidx{:}) = line([X-halfboxwidth X+halfboxwidth],[md md]);
            end
            set(obj.handles.medianLines(subidx{:}),'linestyle','-');
            
            if obj.sampleSize % if sample size requested
                [~, halfboxwidth] = obj.calcBoxWidths(subidx);
                % set x and y offsets
                xoffset = 0.15*halfboxwidth;
                yoffset = (xoffset./(max(obj.xticks)-min(obj.xticks))).*(max(obj.y(:))-min(obj.y(:)));
                gOffset = obj.getOffset(subidx);
                % make text
                obj.handles.samplesTxt(subidx{:}) = text(obj.xticks(subidx{2})+gOffset+halfboxwidth-xoffset,...
                    obj.statistics.PU(subidx{:})-yoffset,...
                    num2str(sum(~isnan(obj.y(subidxAll{:})))),...
                    'horizontalalignment','right','verticalalignment','top');
            end
            
        end
        
        function drawWhiskers(obj,subidx)
        %DRAWWHISKERS draw the whiskers
            
            try % to delete the handles first
                delete(obj.handles.upperWhiskers(subidx{:}));
            catch
            end
            try
                delete(obj.handles.lowerWhiskers(subidx{:}));
            catch
            end
            try
                delete(obj.handles.upperWhiskerTips(subidx{:}));
            catch
            end
            try
                delete(obj.handles.lowerWhiskerTips(subidx{:}));
            catch
            end
            
            % vertices
            X = obj.xticks(subidx{2}) + obj.getOffset(subidx);
            [~, halfboxwidth] = calcBoxWidths(obj,subidx);
            
            % LQ
            obj.handles.lowerWhiskers(subidx{:}) = ...
                line([X X],[obj.statistics.min(subidx{:}) obj.statistics.PL(subidx{:})],'Parent',obj.handles.axes);
            % UQ
            obj.handles.upperWhiskers(subidx{:}) = ...
                line([X X],[obj.statistics.max(subidx{:}) obj.statistics.PU(subidx{:})],'Parent',obj.handles.axes);
            if ~obj.showViolin
                % whisker tips
                obj.handles.lowerWhiskerTips(subidx{:}) = ...
                    line([X-0.5*halfboxwidth X+0.5*halfboxwidth],[obj.statistics.min(subidx{:}) obj.statistics.min(subidx{:})],...
                    'linestyle','-','Parent',obj.handles.axes);
                obj.handles.upperWhiskerTips(subidx{:}) = ...
                    line([X-0.5*halfboxwidth X+0.5*halfboxwidth],[obj.statistics.max(subidx{:}) obj.statistics.max(subidx{:})],...
                    'linestyle','-','Parent',obj.handles.axes);
            end
            
        end
        
        function drawOutliers(obj,subidx)
        %DRAWOUTLIERS draw the outliers
            
            try % to delete the handles first
                delete(obj.handles.outliers(subidx{:}));
            catch
            end
            
            if obj.showOutliers
                if ~isempty(obj.statistics.outliers{subidx{:}}) % don't bother if no data
                    % initial plotting vertices
                    [~, halfboxwidth] = obj.calcBoxWidths(subidx);
                    X = obj.xticks(subidx{2}) + obj.getOffset(subidx);
                    % add a random x offset based on data distribution
                    % use full data to calculate offset
                    subidxAll = subidx;
                    subidxAll{1} = ':';
                    yScatter = obj.y(subidxAll{':'});
                    randOffset = obj.xOffset(yScatter);
                    ix = obj.statistics.outliers_IX(subidxAll{:});
                    ix = ix(~isnan(yScatter));
                    randOffset = randOffset(ix);
                    xScatter = X + (0.8.*halfboxwidth.*randOffset);
                    obj.handles.outliers(subidx{:}) = scatter(xScatter,obj.statistics.outliers{subidx{:}},...
                        'Parent',obj.handles.axes);
                end
            end
            
        end
        
        function drawScatter(obj,subidx,gidx,subidxAll)
        %DRAWSCATTER draw the scatter overlay
            
            try % to delete the handles first
                delete(obj.handles.scatters(subidx{:}));
            catch
            end
            
            if obj.showScatter
                
                % get non-outliers
                IX = ~obj.statistics.outliers_IX(subidxAll{:});
                subidx4 = [{IX} subidx(2:end)];
                
                % get data
                X = obj.xticks(subidx{2}) + obj.getOffset(subidx);
                Y = obj.y(subidx4{:});
                
                % calculate x offsets
                [~,halfboxwidth] = calcBoxWidths(obj,subidx);
                if obj.showViolin
                    halfboxwidth = max(halfboxwidth, obj.calcHalfViolinWidth(subidx));
                end
                xScatter = X + (0.8.*halfboxwidth.*obj.xOffset(Y));
                % plot
                yScatter = Y(~isnan(Y));
                obj.handles.scatters(subidx{:}) = scatter(xScatter,yScatter,obj.scatterMarker{gidx{:}});
                set(obj.handles.scatters(subidx{:}),'Parent',obj.handles.axes);
            end
            
        end
        
        function drawMean(obj,subidx)
        %DRAWMEAN drawn the means
            
            try % to delete the handles first
                delete(obj.handles.means(subidx{:}));
            catch
            end
            
            if obj.showMean
                % x position
                X = obj.xticks(subidx{2}) + obj.getOffset(subidx);
                % plot
                obj.handles.means(subidx{:}) = plot(X,obj.statistics.mean(subidx{:}),'Parent',obj.handles.axes);
            end
            
        end
        
        function drawAddPrctiles(obj,subidx)
        %DRAWADDPRCTILES draw additional percentiles
            
            try
                delete(obj.handles.addPrctiles(subidx{:}))
            catch
            end
            
            if ~isempty(obj.addPrctiles)
                subidx2 = subidx;
                % x position
                X = obj.xticks(subidx{2}) + obj.getOffset(subidx);
                % plot
                h = zeros(length(obj.addPrctiles),1);
                for m = 1:length(obj.addPrctiles)
                    subidx2{1} = m;
                    Y = obj.statistics.addPrctiles(subidx2{:});
                    h(m) = plot(X,Y);
                    text(X,Y,obj.addPrctilesLabels{m},...
                        'horizontalAlignment','left',...
                        'verticalAlignment','bottom')
                end
                subidxAll = subidx;
                subidxAll{1} = ':';
                obj.handles.addPrctiles(subidxAll{:}) = h;
            end
            
        end
        
        function drawGroupGraphics(obj,subidx,gidx)
        %DRAWGROUPGRAPHICS set graphics options for each group
        
            % violin colors
            if obj.showViolin
                set(obj.handles.violin(subidx{:}),...
                    'FaceColor',obj.violinColor{gidx{:}},...
                    'EdgeColor',obj.lineColor{gidx{:}});
            end
        
            % box colors
            set(obj.handles.box(subidx{:}),...
                'FaceColor',obj.boxColor{gidx{:}},...
                'EdgeColor',obj.lineColor{gidx{:}});
            
            % notch line
            if obj.notchLine
                set(obj.handles.notchLowerLines(subidx{:}),...
                    'linestyle',obj.notchLineStyle{gidx{:}},'color',obj.notchLineColor{gidx{:}});
                set(obj.handles.notchUpperLines(subidx{:}),...
                    'linestyle',obj.notchLineStyle{gidx{:}},'color',obj.notchLineColor{gidx{:}});
            end
            
            % median line
            set(obj.handles.medianLines(subidx{:}),'color',obj.medianColor{gidx{:}});
            
            % whiskers
            set([obj.handles.lowerWhiskers(subidx{:}) obj.handles.upperWhiskers(subidx{:})],...
                'color',obj.lineColor{gidx{:}});
            set([obj.handles.lowerWhiskers(subidx{:}) obj.handles.upperWhiskers(subidx{:})],...
                'linestyle',obj.lineStyle{gidx{:}});
            if ~obj.showViolin
                set([obj.handles.lowerWhiskerTips(subidx{:}) obj.handles.upperWhiskerTips(subidx{:})],...
                    'color',obj.lineColor{gidx{:}});
            end
            
            % outliers
            if ~isempty(obj.statistics.outliers{subidx{:}}) && obj.showOutliers % don't bother if no data
                set(obj.handles.outliers(subidx{:}),...
                    'marker',obj.symbolMarker{gidx{:}},'MarkerEdgeColor',obj.symbolColor{gidx{:}});
            end
            
            % scatter
            if obj.showScatter
                set(obj.handles.scatters(subidx{:}),'MarkerEdgeColor',obj.scatterColor{gidx{:}});
            end
            
            % means
            if obj.showMean
                set(obj.handles.means(subidx{:}),...
                    'color',obj.meanColor{gidx{:}},...
                    'marker',obj.meanMarker{gidx{:}});
            end
            
            % additional percentiles
            if ~isempty(obj.addPrctiles)
                subidx2 = subidx;
                for m = 1:length(obj.addPrctiles)
                    subidx2{1} = m;
                    set(obj.handles.addPrctiles(subidx2{:}),...
                        'color',obj.addPrctilesColors{m},...
                        'marker',obj.addPrctilesMarkers{m});
                end
            end
            
        end
        
        function drawGlobalGraphics(obj)
        %DRAWGLOBALGRAPHICS set global graphics options
        
            % violin
            if obj.showViolin
                set(obj.handles.violin, ...
                    'LineWidth',obj.lineWidth, ...
                    'FaceAlpha',obj.violinAlpha);
            end
            
            % box colors
            set(obj.handles.box,...
                'FaceAlpha',obj.boxAlpha,...
                'LineWidth',obj.lineWidth);
            
            % sample size
            if obj.sampleSize
                set(obj.handles.samplesTxt,'fontsize',obj.sampleFontSize);
            end
            
            % notch line
            if obj.notchLine
                set(obj.handles.notchLowerLines,'linewidth',obj.lineWidth);
                set(obj.handles.notchUpperLines,'linewidth',obj.lineWidth);
            end
            
            % median line
            set(obj.handles.medianLines,'linewidth',obj.lineWidth);
            
            % whiskers
            set([obj.handles.lowerWhiskers obj.handles.upperWhiskers],...
                'linewidth',obj.lineWidth);
            if ~obj.showViolin
                set([obj.handles.lowerWhiskerTips obj.handles.upperWhiskerTips],...
                    'linewidth',obj.lineWidth);
            end
            
            % outliers
            if obj.showOutliers && isfield(obj.handles,'outliers')
                set(findobj(obj.handles.outliers,'Type','Scatter'),...
                    'SizeData',obj.outlierSize);
            end
            
            % scatter
            if obj.showScatter && isfield(obj.handles,'scatters')
                set(obj.handles.scatters,'SizeData',obj.scatterSize);
                if obj.scatterAlpha<1
                    try % older versions do not support MarkerEdgeAlpha
                        set(obj.handles.scatters,'MarkerEdgeAlpha',obj.scatterAlpha);
                    catch
                        warning('This version of Matlab does not support scatter marker transparency.')
                    end
                end
            end
            
            % means
            if obj.showMean && isfield(obj.handles,'means')
                set(obj.handles.means,'markersize',obj.meanSize);
            end
            
            % set group label font size
            if isfield(obj.handles,'groupsTxt')
                set(obj.handles.groupsTxt,'FontSize',obj.groupLabelFontSize);
            end
            
            % set additional percentile font size
            if ~isempty(obj.addPrctiles) && isfield(obj.handles,'addPrctilesTxt')
                set(obj.handles.addPrctilesTxt,'FontSize',obj.addPrctilesTxtSize);
            end
            
        end
        
        function drawXsepatator(obj)
        %DRAWXSEPARATOR draw the x separator lines
        
            try % to delete the handles first
                delete(obj.handles.xseps);
            catch
            end
            
            if obj.xSeparator
                xlines = obj.xticks(1:end-1)+0.5.*diff(obj.xticks); % line positions
                obj.handles.xseps = zeros(size(xlines)); % handles
                for n = 1:length(xlines) % draw lines
                    obj.handles.xseps(n) = line([xlines(n) xlines(n)],get(obj.handles.axes,'ylim'),'color',.5+zeros(1,3));
                end
            end
        end
        
        function drawStyle(obj)
        %DRAWSTYLE draw the style options
            
            try % to delete the handles first
                delete(obj.handles.groupsTxt);
            catch
            end
            
            % get outermost box widths (correct xlim for this)
            subidxFirst = cell(1,length(obj.outDims));
            [subidxFirst{:}] = ind2sub(obj.outDims, 1);
            [~,halfboxwidthFirst] = calcBoxWidths(obj,subidxFirst);
            subidxLast = cell(1,length(obj.outDims));
            [subidxLast{:}] = ind2sub(obj.outDims, prod(obj.outDims));
            [~,halfboxwidthLast] = calcBoxWidths(obj,subidxLast);
            
            % set default x axis
            set(obj.handles.axes,'xtick',obj.xticks,...
                'xlim',[min(obj.xticks)-0.5*obj.diffx-halfboxwidthFirst max(obj.xticks)+0.5*obj.diffx+halfboxwidthLast],...
                'xticklabel',obj.xlabels);
            
            switch lower(obj.style)
                case 'normal'
                    % do nothing
                case 'hierarchy'

                    % check grouplabels options
                    if ~isempty(obj.groupLabels)
                        assert(isvector(obj.groupLabels) && iscell(obj.groupLabels),...
                            'iosr:boxPlot:groupLabelsType', ...
                            'The GROUPLABELS option should be a cell vector');
                        if ~isempty(obj.outDims(3:end))
                            assert(isequal(obj.outDims(3:end),cellfun(@length,obj.groupLabels)),...
                                'iosr:boxPlot:groupLabelsSize', ...
                                ['The GROUPLABELS option should be a cell vector; ' ...
                                'the Nth element should contain a vector of length SIZE(Y,N+2)'])
                        else
                            assert(cellfun(@length,obj.groupLabels)==1,...
                                'iosr:boxPlot:groupLabelsSize', ...
                                'The GROUPLABELS option should have length 1 if no grouping dimensions are specified.')
                        end
                    end

                    % number of label rows
                    labelrows = length(obj.outDims)-1;

                    ylim = get(obj.handles.axes,'ylim'); % need y limitis

                    if isnumeric(obj.groupLabelHeight)
                        labelheight = obj.groupLabelHeight;
                    else
                        % determine label height
                        labelheight = obj.calcLabelHeight(); % y range in units calculated based on normalised height
                    end

                    % determine y centre points of labels based on how many labels
                    % and their height
                    ymids = ylim(1) - (labelrows:-1:1).*labelheight + 0.5.*labelheight;

                    % clear axis tick labels
                    set(obj.handles.axes,'xticklabels',[]);

                    % put labels on plot
                    totalticks = prod(obj.outDims(2:end)); % total number of ticks across
                    allticklocs = sort(obj.groupXticks(:)); % their locations
                    obj.handles.groupsTxt = [];
                    for c = labelrows:-1:1 % work down hierarchy
                        Y = ymids(c); % y location
                        nskip = totalticks/prod(obj.outDims(2:c+1)); % skip through tick locations based on hierarchy
                        first = 1:nskip:totalticks; % first tick in group
                        last = nskip:nskip:totalticks; % last tick in group
                        labellocs = allticklocs(first) + 0.5.*(allticklocs(last)-allticklocs(first)); % centre point for label location
                        % work through each label location
                        for n = 1:length(labellocs)
                            if c==1 % special case: x axis ticks
                                gLabel = obj.xlabels{n};
                            else % every other group
                                if isempty(obj.groupLabels) % auto increment label
                                    gLabel = num2str(mod(n-1,obj.outDims(c+1))+1);
                                else % get from input
                                    gLabel = obj.groupLabels{c-1}(mod(n-1,obj.outDims(c+1))+1); % cyclic index into labels
                                    if isnumeric(gLabel) % convert numbers to strings
                                        gLabel = num2str(gLabel);
                                    end
                                end
                            end
                            % place actual label
                            t = text(labellocs(n),Y,gLabel,'HorizontalAlignment','center','VerticalAlignment','middle');
                            % store info about placement
                            setappdata(t,'rows',labelrows);
                            setappdata(t,'row',c);
                            % add to handles
                            obj.handles.groupsTxt = [obj.handles.groupsTxt t];
                        end

                    end
                otherwise
                    error('iosr:boxPlot:unknownStyle','Unkown ''style'' option.')
            end
            
        end
        
        function drawLegend(obj)
        %DRAWLEGEND draw the legend
        
            try % to delete the handles first
                delete(obj.handles.legend);
            catch
            end
            
            if obj.showLegend
                % dimensions of group data
                gDims = cellfun(@length,obj.groupLabels);
                if isempty(gDims)
                    gDims = 1;
                end
                subidx = cell(1,length(gDims));
                % going to change order of handles
                orderidx = cell(1,length(gDims)+1);
                order = zeros(prod(gDims),1);
                % preallocate legend labels
                lgndstr = cell(prod(gDims),1);
                for n = 1:prod(gDims)
                    % use linear index to get group data
                    [subidx{:}] = ind2sub([gDims 1], n);
                    % order looks at trailing dimensions first
                    [orderidx{:}] = ind2sub(fliplr([gDims 1]), n);
                    orderidx = fliplr(orderidx);
                    order(n) = sub2ind([gDims 1],orderidx{:});
                    % create labels
                    start = true; % create start of string label
                    for m = length(subidx):-1:1 % work from end through dimensions
                        str = obj.ensureString(obj.groupLabels{m}(subidx{m}));
                        if start % initialise label
                            lgndstr{n} = str;
                            start = false;
                        else % prepend string label
                            lgndstr{n} = [str ', ' lgndstr{n}];
                        end
                    end
                end
                % choose target for the legend
                if ~obj.arrayElementsEqual(obj.boxColor)
                    legendTarget = obj.handles.box;
                    target = 'box';
                elseif ~obj.arrayElementsEqual(obj.medianColor)
                    legendTarget = obj.handles.medianLines;
                    target = 'medianLines';
                elseif obj.showViolin && ~obj.arrayElementsEqual(obj.violinColor)
                    legendTarget = obj.handles.violin;
                    target = 'violin';
                else
                    error('iosr:boxPlot:legend','The legend uses the box, median line colors, or violins to create legend entries. However, these items appear to be identical, which would make the legend impossible to interpret. Change the ''boxColor'', ''medianColor'', or ''violinColor'' option.');
                end
                % create legend in specified order
                if ~isempty(legendTarget)
                    [obj.handles.legend, icons] = legend(legendTarget(1,1,order),lgndstr(order));
                    set(obj.handles.legend,'location','best');
                    % add alpha to legend entries
                    PatchInLegend = findobj([obj.handles.legend; icons(:)], 'type', 'patch');
                    if ~isempty(PatchInLegend)
                        switch lower(target)
                            case 'box'
                                set(PatchInLegend, 'facealpha', obj.boxAlpha);
                            case 'violin'
                                set(PatchInLegend, 'facealpha', obj.violinAlpha);
                        end
                    end
                end
            end
        end
        
        function createTheme(obj)
        %CREATETHEME create the color/line theme
            
            switch lower(obj.theme) % choose scheme
                case 'colorall'
                    boxalpha = 0.4;
                    boxcolor = obj.themeColors;
                    linecolor = 'k';
                    linestyle = '-';
                    meancolor = 'k';
                    mediancolor = 'k';
                    notchlinecolor = obj.themeColors;
                    notchlinestyle = ':';
                    scattercolor = obj.themeColors;
                    symbolcolor = obj.themeColors;
                case {'colorlines', 'graylines'}
                    if strcmpi(obj.theme,'graylines')
                        warning('iosr:boxPlot:graylines', 'The ''graylines'' theme is deprecated. Use the ''colorlines'' theme instead, and specify ''themeColors'' as @gray.')
                    end
                    boxalpha = obj.boxAlpha;
                    boxcolor = 'none';
                    linecolor = 'k';
                    linestyle = '-';
                    meancolor = 'k';
                    mediancolor = obj.themeColors;
                    notchlinecolor = obj.themeColors;
                    notchlinestyle = ':';
                    scattercolor = obj.themeColors;
                    symbolcolor = obj.themeColors;
                case {'colorboxes', 'grayboxes'}
                    if strcmpi(obj.theme,'grayboxes')
                        warning('iosr:boxPlot:grayboxes', 'The ''grayboxes'' theme is deprecated. Use the ''colorboxes'' theme instead, and specify ''themeColors'' as @gray.')
                    end
                    boxalpha = obj.boxAlpha;
                    boxcolor = obj.themeColors;
                    linecolor = 'k';
                    linestyle = '-';
                    meancolor = 'k';
                    mediancolor = 'k';
                    notchlinecolor = 'k';
                    notchlinestyle = ':';
                    scattercolor = [.5 .5 .5];
                    symbolcolor = 'k';
                otherwise
                    if ~strcmpi(obj.theme,'default')
                        warning(['Unknown theme ''' obj.theme ''' specified. Using default'])
                    end
                    boxalpha = 1;
                    boxcolor = 'none';
                    linecolor = 'k';
                    linestyle = '-';
                    meancolor = obj.themeColors;
                    mediancolor = obj.themeColors;
                    notchlinecolor = 'k';
                    notchlinestyle = ':';
                    scattercolor = [.5 .5 .5];
                    symbolcolor = obj.themeColors;
            end
            % turn off legend while changing parameters
            legendState = obj.showLegend;
            obj.showLegend = false;
            % set parameters
            obj.boxAlpha = boxalpha;
            obj.boxColor = boxcolor;
            obj.medianColor = mediancolor;
            obj.lineColor = linecolor;
            obj.lineStyle = linestyle;
            obj.meanColor = meancolor;
            obj.notchLineColor = notchlinecolor;
            obj.notchLineStyle = notchlinestyle;
            obj.scatterColor = scattercolor;
            obj.symbolColor = symbolcolor;
            % turn legend back on if it was on
            obj.showLegend = legendState;
        end
        
        function calculateStats(obj)
        %CALCULATESTATS calculate the statistic for the box plot
        
            calculateStats@iosr.statistics.statsPlot(obj);
            
            % calculate stats
            obj.statistics.percentile = obj.percentile; % percentile
            outsize = size(obj.statistics.median);
            obj.statistics.addPrctiles = zeros([length(obj.addPrctiles) outsize(2:end)]);
            obj.statistics.PL = zeros(outsize); % lower percentile
            obj.statistics.PU = zeros(outsize); % upper percentile
            subidx = cell(1,length(obj.outDims));
            for n = 1:prod(obj.outDims)
                [subidx{:}] = ind2sub(obj.outDims, n);
                subidxAll = subidx;
                subidxAll{1} = ':';
                if ~isempty(obj.addPrctiles)
                    obj.statistics.addPrctiles(subidxAll{:}) = iosr.statistics.quantile(obj.y(subidxAll{:}),obj.addPrctiles(:)./100,[],obj.method,obj.weights(subidxAll{:}));
                end
                obj.statistics.PL(subidx{:}) = iosr.statistics.quantile(obj.y(subidxAll{:}),min(obj.percentile)/100,[],obj.method,obj.weights(subidxAll{:}));
                obj.statistics.PU(subidx{:}) = iosr.statistics.quantile(obj.y(subidxAll{:}),max(obj.percentile)/100,[],obj.method,obj.weights(subidxAll{:}));
            end

            % check for notches extending beyond box
            if (any(obj.statistics.notch_u(:)>obj.statistics.Q3(:)) || any(obj.statistics.notch_l(:)<obj.statistics.Q1(:))) && (obj.notch || obj.notchLine)
                warning('Notch extends beyond quartile. Try setting ''notch'' or ''notchLine'' to false')
            end
            
        end
        
        function setXprops(obj)
        %SETXPROPS set x axis properties
            
            if isnumeric(obj.x) % x is numeric
                obj.xlabels = strtrim(cellstr(num2str(obj.x(:))));
                switch lower(obj.xSpacing) % choose spacing
                    case 'x'
                        obj.xticks = obj.x;
                        obj.diffx = min(diff(obj.x));
                        if isempty(obj.diffx)
                            obj.diffx = 1;
                        end
                    case 'equal'
                        obj.xticks = 1:length(obj.x);
                        obj.diffx = 1;
                    otherwise
                        error('iosr:boxPlot:xSpacing','Unknown xSpacing option specified.')
                end
            else % x is not numeric
                obj.xlabels = obj.x;
                obj.xticks = 1:length(obj.x);
                obj.diffx = 1;
            end
            obj.groupRange = obj.groupWidth*obj.diffx;
        end
        
        function out = groupOption(obj,option,name)
        %GROUPOPTION convert parameter to format for plotting

            gDims = obj.groupDims;
            if length(gDims)==1
                gDims = [gDims 1];
            end

            % checks
            if ischar(option)
                option = cellstr(option); % ensure string is cell array
            elseif isnumeric(option)
                option = {option};
            end
            if isvector(option)
                option = option(:);
            end

            % pre-allocate output
            out = cell(gDims);

            % create output
            if length(option)==1 % repeat if only one specified
                for n = 1:prod(gDims)
                    out(n) = option;
                end
            elseif isequal(gDims,size(option)) % put into cell array
                for n = 1:prod(gDims)
                    assert(ischar(option{n}) || size(option{n},2)==3, 'iosr:boxPlot:groupOptionFormat', ['Option ''' name ''' is not in the correct format.'])
                    out(n) = option(n);
                end
            else
                error('iosr:boxPlot:groupOptionSize','Option ''%s'' is not the correct size.',name)
            end

        end

        function out = checkColor(obj,color,usage)
        %checkColor convert 'auto' color option to RGB

            gDims = obj.groupDims;
            out = color; % pass input to output
            
            if ischar(color)
                if strcmpi(color,'auto')
                    color = @lines;
                end
            end
            
            if isa(color,'function_handle')
                obj.checkColorHandle(color);
                out = cell([gDims 1]);
                if any(strcmpi(char(color),{'hot','gray','bone','copper','pink'}))
                    N = 512;
                    cmapTest = color(N);
                    YL = 0.2989.*cmapTest(:,1) + 0.5870*cmapTest(:,2) + 0.1140.*cmapTest(:,3); % luminance
                    switch lower(usage)
                        case 'fill'
                            % ensure minimum luminance of 0.33
                            firstLumaIX = find(YL>0.33,1,'first');
                            if numel(out)<=N
                                cmap = cmapTest(round(linspace(firstLumaIX,N,numel(out))),:); % sample color map
                            else
                                pad = ceil((firstLumaIX/N)*numel(out)); % pad for the black region
                                cmap = color(numel(out)+pad); % colormap
                                cmap = cmap(pad+1:end,:); % remove padded region
                            end
                        case 'line'
                            % ensure maximum luminance of 0.66
                            lastLumaIX = find(YL<0.66,1,'last');
                            if numel(out)<=N
                                cmap = cmapTest(round(linspace(1,lastLumaIX,numel(out))),:); % sample color map
                            else
                                pad = ceil(((N-lastLumaIX)/N)*numel(out)); % pad for the black region
                                cmap = color(numel(out)+pad); % colormap
                                cmap = cmap(1:end-pad,:); % remove padded region
                            end
                    end
                else
                    cmap = color(numel(out)); % colormap
                end
                for n = 1:numel(out)
                    out{n} = cmap(n,:);
                end
            end

        end
        
        function [boxwidth, halfboxwidth] = calcBoxWidths(obj,subidx)
        %CALCBOXWIDTHS calculate the box width
            
            % size of boxes
            if ischar(obj.boxWidth)
                boxwidth =  0.75*(obj.groupRange/prod(obj.groupDims));
            else
                boxwidth = obj.boxWidth;
            end
            
            if obj.scaleWidth
                % scale box width according to sample size
                maxN = max(obj.statistics.N(:));
                halfboxwidth = 0.5 * (sqrt(obj.statistics.N(subidx{:})/maxN)*boxwidth);
            else
                % normal box width
                halfboxwidth = boxwidth/2;
            end
            
        end
        
        function halfviolinwidth = calcHalfViolinWidth(obj,subidx)
        %HALFVIOLINWIDTH calculate the violin width
            
            % size of boxes
            if ischar(obj.violinWidth)
                halfviolinwidth =  0.75*(obj.groupRange/prod(obj.groupDims));
            else
                halfviolinwidth = obj.violinWidth;
            end
            
            if obj.scaleWidth
                % scale box width according to sample size
                maxN = max(obj.statistics.N(:));
                halfviolinwidth = 0.5 * (sqrt(obj.statistics.N(subidx{:})/maxN)*halfviolinwidth);
            else
                halfviolinwidth = halfviolinwidth / 2;
            end
            
        end
        
        function labelheight = calcLabelHeight(obj)
        %CALCLABELHEIGHT calculate the label height for the axes

            ylim = get(obj.handles.axes,'ylim');
            plotheight = get(obj.handles.axes,'position'); % plot position
            plotheight = plotheight(4); % plot height in normalised units
            plotrange = abs(diff(ylim)); % y range
            labelheight = plotheight*plotrange*0.05; % y range in units calculated based on normalised height

        end
        
        function offset = getOffset(obj,subidx)
        %GET_OFFSET calculate offset for hierarchical data

            dims = obj.outDims(3:end);
            subidx = subidx(3:end);
            if isempty(dims) || prod(dims)==1
                offset = 0;
            else
                boxspacing = obj.groupRange/prod(obj.groupDims);
                offsetn = (0:prod(dims)-1)-((prod(dims)-1)./2);
                n = sub2ind([dims(end:-1:1) 1],subidx{end:-1:1});
                offset = offsetn(n).*boxspacing;
            end

        end
        
        function val = checkAddPrctileProp(obj,init,default)
        %CHECKADDPRCTILEPROP check and set default values of additional percentile properties
            
            val = init(:);
            currlength = length(val);
            if currlength<length(obj.addPrctiles)
                val = [val; cell(length(obj.addPrctiles)-currlength,1)];
                for m = currlength+1:length(obj.addPrctiles)
                    val(m) = default;
                end
            end
        end
        
        function eventHandler(obj,hProp,eventData)
        %EVENTHANDLER handle axes property change events

            propName = hProp.Name;
            propValue = eventData.AffectedObject.(propName);

            % update separator lines
            if strcmpi(propName,'YLim')
                % correct x separator limits
                if isfield(obj.handles,'xseps')
                    ydata = propValue;
                    if any(isinf(propValue))
                        % inf uses actual data to determine limits
                        ylim = [NaN NaN];
                        children = get(obj.handles.axes,'children');
                        for n = 1:length(children)
                            try
                                d = get(children(n),'YData');
                                ylim(1) = min([ylim(1); d]);
                                ylim(2) = max([ylim(2); d]);
                            catch
                                % ignore objects that don't have YData
                            end
                        end
                        % replace occurences of inf
                        ydata(isinf(propValue)) = ylim(isinf(propValue));
                    end
                    set(obj.handles.xseps,'YData',ydata);
                end
                % move group labels up or down
                if isfield(obj.handles,'groupsTxt')
                    labelheight = obj.calcLabelHeight();
                    for n = 1:length(obj.handles.groupsTxt)
                        % current position
                        pos = get(obj.handles.groupsTxt(n),'Position');
                        % new position
                        rows = getappdata(obj.handles.groupsTxt(n),'rows');
                        row = getappdata(obj.handles.groupsTxt(n),'row');
                        pos(2) = propValue(1) - (rows-row+.5)*labelheight;
                        % set it
                        set(obj.handles.groupsTxt(n),'Position',pos);
                    end
                end
            end

        end
        
    end
    
    methods (Static, Access = private)
        
        function offset = xOffset(y)
        %XOFFSET create an x offset based on data distribution

            y = y(~isnan(y));
            [d, xd] = iosr.statistics.kernelDensity(y);
            d = d./max(d);
            
            try
                maxDisplacement = interp1(xd, d, y);
            catch
                maxDisplacement = zeros(size(y));
            end
            randOffset = rand(size(y));%randperm(numel(y))-1;
            randOffset = (2*randOffset) - 1;
            % randOffset = (2*randOffset./max(randOffset)) - 1;
            offset = randOffset(:) .* maxDisplacement;
            
        end
        
        function str = ensureString(val)
        %ENSURESTRING ensure input value is a string or convert to a string
            
            if isnumeric(val)
                str = num2str(val);
            elseif ischar(val)
                str = val;
            elseif iscellstr(val)
                str = char(val);
            else
                error('iosr:boxPlot:ensureString','Unable to convert data type to string.')
            end
            
        end
        
        function equal = arrayElementsEqual(v)
        %ARRAYELEMENTSEQUAL check if elements in array are the same
            
            x = v(:);
            if iscell(x) && ~iscellstr(x)
                x = cell2mat(x); % convert numeric cell arrays to numeric arrays
            end
            equal = numel(unique(x))==1;
            
        end
        
        function checkColorHandle(fHandle)
        %CHECKCOLORHANDLE check whether function handle reurns valid color output
            
            for N = [1 4 8 16]
                assert(isequal(size(fHandle(N)),[N 3]), 'iosr:boxPlot:colorFhandle', ['The color function ''' char(fHandle) ''' does not appear to return a color array of the correct size.'])
                cOutput = fHandle(N);
                assert(all(cOutput(:)>=0) && all(cOutput(:)<=1), 'iosr:boxPlot:colorFhandleRange', ['The color ''' char(fHandle) ''' function does not appear to return RGB values in the interval [0,1].'])
            end
        end
        
        function option = replaceAuto(option, default, errorMsg)
        %REPLACEAUTO check for and replace auto parameters
            if ischar(option)
                if strcmpi(option, 'auto')
                    option = default;
                else
                    error('iosr:boxPlot:invalidOption',errorMsg)
                end
            end
        end
        
    end
    
end
