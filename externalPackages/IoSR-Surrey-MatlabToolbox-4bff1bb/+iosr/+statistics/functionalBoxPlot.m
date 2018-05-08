classdef (CaseInsensitiveProperties = true) functionalBoxPlot < ...
        iosr.statistics.functionalPlot
%FUNCTIONALBOXPLOT Draw a functional boxplot
% 
%   Use this class to plot a functional boxplot [1]. A functional boxplot
%   considers each function as a single observation. The order statistics
%   are calculated from the band depth (BD) or modified band depth (MBD) of
%   each function; the appearance of the functional boxplot is derived
%   analogously to the traditional boxplot. The plot shows the median
%   function, the 50% central region, and the functional limits
%   ([Q1-1.5*IQR, Q3+1.5*IQR], where Q1 and Q3 are the 25th and 75th
%   percentiles respectively, and IQR is the interquartile range). Outliers
%   (functions exceeding this limit) may also be plotted.
% 
%   FUNCTIONALBOXPLOT operates on the array Y, which is an N-by-M-by-P
%   array, where N functions are stored across the rows of Y, there are M
%   X-values, and there are P functional plots.
%   
%   IOSR.STATISTICS.FUNCTIONALBOXPLOT properties:
%       method                  - Specifies the method used to calculate
%                                 the band depth.
%       outlierLineColor        - Specifies the color of the outlier lines.
%       outlierLineStyle        - Specifies the width of the outlier marker
%                                 line. See 'mainLineWidth' for details of
%                                 how to specify widths. The default is 1.
% 
%   Additional properties are inherited from
%   IOSR.STATISTICS.FUNCTIONALPLOT.
% 
%   These properties can be referenced using dot notation - e.g. H.METHOD
%   where H is an instance of the FUNCTIONALBOXPLOT object - or using
%   the SET and GET methods - e.g. GET(H,'METHOD'). Both methods are
%   case-insensitive.
% 
%   Note that some handles will be empty unless the associated option is
%   turned on.
% 
%   Read-only properties:
%       x                   - The x data.
%       y                   - The y data.
%       statistics          - Structure containing the statistics used
%                             for the box plot. With the exception of
%                             'outliers', each field contains a 1-by-M-by-P
%                             numeric array of values (identical to that
%                             returned by IOSR.STATISTICS.QUANTILE). The
%                             fields are:
%                                 'inner_l'       : the 25th percentile
%                                 'inner_u'       : the 75th percentile
%                                 'main'          : the values of the
%                                                   central line
%                                 'outer_l'       : the lower limits
%                                 'outer_u'       : the upper limits
%                                 'outliers'      : a 1-by-P cell array of
%                                                   outlier functions
% 
%   IOSR.STATISTICS.FUNCTIONALBOXPLOT methods:
%       FUNCTIONALBOXPLOT    - Create the plot.
% 
%   References
% 
%   [1] Sun, Y.; Genton, M. G. (2011). "Functional boxplots". Journal of
%       Computational and Graphical Statistics. 20: 316?334.
%       doi:10.1198/jcgs.2011.09224
% 
%   See also IOSR.STATISTICS.BOXPLOT, IOSR.STATISTICS.QUANTILE,
%   IOSR.STATISTICS.FUNCTIONALPLOT, IOSR.STATISTICS.FUNCTIONALSPREADPLOT.

% Adapted from code originally written by:
%  - Ying Sun: sunwards@stat.tamu.edu
%  - Marc G. Genton: genton@stat.tamu.edu

    properties (AbortSet)
        method = 'mbd'                  % The method used to calculate the quantiles.
        outlierLineColor = 'auto'       % Color of the outlier lines
        outlierLineStyle = '-'          % Width of the outlier lines.
    end
    
    properties (Access = private)
        depths
    end
    
    methods
        
        % constructor
        function obj = functionalBoxPlot(varargin)
        % FUNCTIONALBOXPLOT Draw a functional boxplot.
        %
        %   IOSR.STATISTICS.FUNCTIONALBOXPLOT(Y) produces a functional
        %   boxplot of the data in Y. Y should be an N-by-M or N-by-M-by-P
        %   array, where N functions are stored in the rows of Y, there are
        %   M X-values, and there are P functional plots. The plot shows
        %   the median function as the main line, the interquartile range
        %   as the shaded region, lines showing the limits of the data, and
        %   lines showing any outlier functions.
        %
        %   Tabular data can be arranged into the appropriate format using
        %   the IOSR.STATISTICS.TAB2BOX function.
        %
        %   IOSR.STATISTICS.FUNCTIONALBOXPLOT(X,Y) specifies the x-axis
        %   values for the plot. X should be an M-length vector. The
        %   default is 1:M.
        %
        %   IOSR.STATISTICS.FUNCTIONALBOXPLOT(...,'PARAMETER',VALUE)
        %   allows the plotting options to be specified when the plot is
        %   constructed.
        %
        %   IOSR.STATISTICS.FUNCTIONALBOXPLOT(AX,...) creates the plot
        %   in the axes specified by AX.
        %   
        %   Example
        %
        %       % generate random data
        %       y = cat(3, randn(100,21), 0.25+randn(100,21));
        %       x = 0:20;
        %    
        %       % Draw a functional box plot for the first function
        %       figure
        %       iosr.statistics.functionalBoxPlot(x, y(:,:,1));
        %       title('Functional boxplot.')
        %       axis tight
        %       box on
        
            start = obj.getXY(varargin{:});
            
            % check input is valid size
            assert(ndims(obj.y) <= 3 && ndims(obj.y) >= 2, ...
                'iosr:functionalBoxPlot:invalidY', ...
                'Y must be a two- or three-dimensional array.');
            
            % set the properties of the plot
            obj.whiskers = true;
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
        
        function set.method(obj, val)
            obj.method = val;
            obj.calculateStats();
            obj.draw();
        end
        
        function val = get.outlierLineColor(obj)
            val = obj.parseColor(obj.outlierLineColor);
        end
        
        function val = get.outlierLineStyle(obj)
            val = obj.parseProps(obj.outlierLineStyle, false);
        end
        
    end
    
    methods (Access = protected)
        
        % calculate the statistics used in the plot
        function calculateStats(obj)
            
            obj.depths = obj.bandDepth();
            nCurves = size(obj.y, 1);
            
            % Calculate specific stats
            outsize = [1 obj.ydims(2:end)];
            outDimsTemp = [1 1 obj.outDims(3:end)];
            obj.statistics.main = zeros(outsize);
            obj.statistics.inner_u = zeros(outsize);
            obj.statistics.inner_l = zeros(outsize);
            obj.statistics.outer_u = zeros(outsize);
            obj.statistics.outer_l = zeros(outsize);
            subidx = cell(1,length(obj.outDims));
            for n = 1:prod(outDimsTemp)
                [subidx{:}] = ind2sub(outDimsTemp, n);
                subidx{2} = ':';
                subidxAll = subidx;
                subidxAll{1} = ':';
                
                data = obj.y(subidxAll{:});
                
                % main statistics
                [~,index] = sort(obj.depths,'descend');
                m = ceil(nCurves*0.5);
                center = data(index(1:m),:);
                inf = min(center);
                sup = max(center);
                dist = 1.5*(sup-inf);
                upper = sup+dist;
                lower = inf-dist;
                
                % outliers
                outly = sum(or(data'<=lower'*ones(1,nCurves), data'>=upper'*ones(1,nCurves)))';
                outpoint = find(outly);
                outliers = data(outpoint,:);
                good = data;
                good(outpoint,:)=[];
                maxcurve = max(good);
                mincurve = min(good);
                
                % parent class uses these fields for plotting
                obj.statistics.inner_u(subidx{:}) = sup;
                obj.statistics.inner_l(subidx{:}) = inf;
                obj.statistics.outer_u(subidx{:}) = maxcurve;
                obj.statistics.outer_l(subidx{:}) = mincurve;
                obj.statistics.main(subidx{:}) = data(index(1),:);
                obj.statistics.outliers{n} = outliers;
                
            end
            
        end
        
        % Draw the plot
        function draw(obj)
            
            if isfield(obj.handles,'axes')

                axes(obj.handles.axes);
                cla;
                hold on;
                
                % parent class does most of the plotting
                draw@iosr.statistics.functionalPlot(obj);
                
                % plot outlier functions
                if length(obj.outDims) > 2
                    nlines = obj.outDims(3);
                else
                    nlines = 1;
                end
                if obj.showOutliers
                    for n = nlines:-1:1
                        % draw inner patch
                        if ~isempty(obj.statistics.outliers{n})
                            obj.handles.outliers{n} = line(obj.x, obj.statistics.outliers{n}, ...
                                'linestyle', obj.outlierLineStyle{n}, ...
                                'linewidth', obj.outlierLineWidth{n}, ...
                                'color', obj.outlierLineColor{n} ...
                            );
                        end
                    end
                else
                    try
                        delete(obj.handles.outliers)
                    catch
                    end
                end
                
            end

        end
        
    end
    
    methods (Access = private)
        
        % calculate band depth
        function depth = bandDepth(obj)
            
            switch lower(obj.method)
                case 'bd2'
                    depth = obj.BD2();
                case 'bd3'
                    depth = obj.BD3();
                case 'mbd'
                    depth = obj.MBD();
                otherwise
                    error('iosr:functionalBoxPlot:unknownMethod',['Unknown mode ''' obj.method ''''])
            end

        end


        function contg = MBD(obj)
            % calculates the generalized band depth of a set of data
            n = size(obj.y,1); % size of the data matrix
            cont = zeros(n,1,obj.outdims(3));
            for p = 1:obj.outdims(3)
                for i=1:(n-1)
                    for j=(i+1):(n) % consider all possible pairs of functions
                        cont(:,:,p) = cont(:,:,p)+obj.a(obj.y(:,:,p), [i j]);
                    end
                end
            end
            contg=cont/obj.combinat(n,2);
        end

        function contg = BD3(obj)
            % calculate the band depth with J=3.
            n=size(obj.y,1); 
            cont=zeros(n,1,obj.outdims(3));
            % Select three observations from the sample in all the possible ways.
            for p = 1:obj.outdims(3)
                for i=1:(n-2)
                   for j=(i+1):(n-1)
                      for k=(j+1):n
                         cont(:,:,p) = cont(:,:,p)+obj.estaEntre(obj.y(:,:,p), [i j k])';   % In this subfunction we check which observations from the sample is inside the band delimeted by observations i,j and k.           

                      end
                   end
                end	
            end
            contg=cont/obj.combinat(n,3);
        end

        function contg = BD2(obj)
            % calculate the band depth of every observation in the matrix
            n=size(obj.y,1);
            cont=zeros(n, 1, obj.outdims(3));
            for p = 1:obj.outdims(3)
                for i=1:((n+2)/3)
                   for j=(i+1):(n) % choose pairs of indexes in all the possible ways.
                      cont(:,:,p) = cont(:,:,p)+obj.estaEntre(obj.y(:,:,p), [i j])';
                   end
                end
            end
            contg =cont/obj.combinat(n,2);
        end
        
    end
    
    methods (Static, Access = private)
        
        function resultado = a(data, v) 
            n = size(data,1);
            p = size(data,2);
            Z = data;
            inf=(min(Z(v,:)))';
            sup=(max(Z(v,:)))';
            resul=sum((and(Z'<=sup*ones(1,n),Z'>=inf*ones(1,n))));
            % Proportion of coordinates of each observation from the sample
            % that is inside the band delimited by pairs v of functions from the sample.
            resultado=(resul/p)';
        end
        
        function resultados = estaEntre(data, v)
            [n,p]=size(data);                      
            Z=data;
            inf=min(Z(v,:))';
            sup=max(Z(v,:))';
            resultados=sum((and(Z'<=sup*ones(1,n),Z'>=inf*ones(1,n))))==p;
        end
        
        function combinat = combinat(n,p)
            if n<p
            combinat=0;
            else
               combinat = nchoosek(n,p);
            end
        end
        
    end
    
end
